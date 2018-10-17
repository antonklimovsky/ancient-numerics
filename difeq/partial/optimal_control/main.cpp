//----------------------------------------------
// Solves folowing optimal control task
//
// \begin{equation}
// \label{eq:func}
// I(x,u) := \frac{1}{2}\int_0^T{[(C{x(t)-z_0(t)},
//           x(t)-z_0(t))+(Nu(t),u(t))]dt} \\
// \end{equation}    
//
// \begin{eqnarray}
// \label{eq:pde}
// & \dot{x} = Ax+Bu, \ t>0 & \\
// & x(0) = x_0. & \nonumber 
// \end{eqnarray}    
//
// Author: A. Klimovsky, root@ludus.kharkiv.com
//----------------------------------------------
#include <iostream.h>
#include <math.h>

#include "lin_alg.h"
#include "gauss.h"

//-------------------------------------------
// Global data
//-------------------------------------------
//----------- input parameters
double T = M_PI;
double eps = 0.05;
int nn = 2;            // system state vector dimension
int mm = 2;            // control vector dimension
long int nnn = 20;     // output mesh cardinality
long int p;            // current mesh cardinality
double h;              // step between nodes of the mesh
matrix c;              
matrix n;
matrix a;
matrix b;
vector x0;
//----------- calculated variables
vector xOld;
matrix d;
matrix nInv;
matrix aTransp;
matrix* q;
vector* r;
vector* x;
vector* u;
double functional;

//---------------------------------------------------------------
// Runge-Kutta method implementation
//
// @param T* x -- should point to a sufficient memory area BEFORE
//                invocation of the function.
//---------------------------------------------------------------
template<class T>
void rungeKutt(
    double h,
    long int n,
    double t0,
    T x0,
    T (*f)(double, T),
    T* x
)
{
    T k0, k1, k2, k3;
    int i;
    double t;

    if (h > 0) {
        x[0] = x0;
        for (i = 1, t = t0; i < n; i++, t += h) {
            k0 = h*(*f)(t, x[i-1]);
            k1 = h*(*f)(t+h/2, x[i-1]+k0/2);
            k2 = h*(*f)(t+h/2, x[i-1]+k1/2);
            k3 = h*(*f)(t+h, x[i-1]+k2);
            x[i] = x[i-1]+(k0+2*k1+2*k2+k3)/6;
        }
    } else {
        x[n-1] = x0;
        for (i = n-2, t = t0; i >= 0; i--, t += h) {
            k0 = h*f(t, x[i+1]);
            k1 = h*f(t+h/2, x[i+1]+k0/2);
            k2 = h*f(t+h/2, x[i+1]+k1/2);
            k3 = h*f(t+h, x[i+1]+k2);
            x[i] = x[i+1]+(k0+2*k1+2*k2+k3)/6;
        }
    }
}

void init()
{
    double myC[] = {
        1, 0, 
        0, 1
    };
    double myN[] = {
        eps*1, eps*0,
        eps*0, eps*1
    };
    double myA[] = {
        0, 1,
       -1, 0
    };
    double myB[] = {
        0, 1,
        1, 2
    };
    double myX0[] = {1, 1};

    c = matrix(nn, nn, myC);
    n = matrix(mm, mm, myN);
    a = matrix(nn, nn, myA);
    b = matrix(mm, mm, myB);
    x0 = vector(nn, myX0);
    nInv = calculate_inverse_matrix(n);
    d = -b*nInv*b.transpose();
    aTransp = a.transpose();
}

vector z0(double t)
{
    vector temp(nn);

    temp[0] = 0;
    temp[1] = 0;

    return temp;
}

//--------------------
// Utility function f
//--------------------
vector f(double t)
{
    return -c*z0(t);
}

long int fix(double t)
{
    long int result;

    result = t/h;
    if (result < 0)
        result = 0;
    if (result >= p)
        result = p-1;
        
    return result;
}

//----------------------------------
// r.h.s. of ODE for q (q'=fq(t,q))
//----------------------------------
matrix fq(double t, matrix q)
{
    return -(q*(a+d*q)+c+aTransp*q);
}

//----------------------------------
// r.h.s. of ODE for r (r'=fr(t,r))
//----------------------------------
vector fr(double t, vector r)
{
    return -((q[fix(t)]*d+aTransp)*r+f(t));
}

//----------------------------------
// r.h.s. of ODE for x (x'=fx(t,x))
//----------------------------------
vector fx(double t, vector x)
{
    return (a+d*q[fix(t)])*x+d*r[fix(t)];
}

void calcQ()
{
    matrix zero(nn, mm);
    int i, j;

    for (i = 0; i < nn; i++)
        for (j = 0; j < mm; j++)
            zero[i][j] = 0;

    q = new matrix[p+1];
    rungeKutt(-h, p+1, T, zero, fq, q);
}

void calcR()
{
    vector zero(nn);
    int i;

    for (i = 0; i < nn; i++)
        zero[i] = 0;

    r = new vector[p+1];
    rungeKutt(-h, p+1, T, zero, fr, r);
}

void calcX()
{
    x = new vector[p+1];
    rungeKutt(h, p+1, 0.0, x0, fx, x);
}

void calcU()
{
    matrix temp;
    long int i;

    temp = -nInv*b.transpose();
    u = new vector[p+1];
    for (i = 0; i < p+1; i++)
        u[i] = temp*(q[i]*x[i]+r[i]);
}

void calcFunctional()
{
    long int i;
    double t;

    functional = 0;
    for (i = 0, t = 0; i < p; i++, t += h)
        functional += ((c*(x[i]-z0(t)))*(x[i]-z0(t))+
                       (n*u[i])*(u[i]))*h;
    functional /= 2;
}

void iterate()
{
    h = T/p;
    calcQ();
    calcR();
    calcX();
}

void output()
{
    long int i;
    long int step;

/*    cout.setf(ios_base::fixed, ios_base::floatfield);*/
    cout.precision(3);

    cout << "Mesh cardinality: " << p << endl;
    step = p/nnn;
    for (i = 0; i < p+1; i += step) {
        cout << "x(" << i*h << ")" << ": " << x[i];
        cout << "u(" << i*h << ")" << ": " << u[i] << endl;
    }
    cout << "Optimal value of the functional: " << functional  << endl;
}

void destroy()
{
    delete [] q;
    delete [] r;
    delete [] x;
    delete [] u;
}

void main()
{
    init();
    p = nnn*10;
    iterate();
    xOld = x[p-1];
    while (1) {
        destroy();
        p <<= 1;
        iterate();
        if (abs(x[p-1]-xOld) > eps)
            xOld = x[p-1];
        else {
            calcU();
            calcFunctional();
            output();
            destroy();
            break;
        }
    }
}