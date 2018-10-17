//-----------------------------------------------------------//
//  Method of Galerkin for solving the boundary problem for  //
//  ordinary differential equations of the 2nd degree.       //
//                                                           //
//  Author: A. Klimovsky                                     //
//-----------------------------------------------------------//
#include <iostream.h>
#include <math.h>

#include "difeq.h"
#include "progonka.h"
#include "integrator.h"
#include "galerkin.h"

Galerkin::Galerkin(DifEq& myEq):
    DifEqSolver(myEq)
{};

inline double Galerkin::node(long int i) // h should be already assigned!
{
    return eq.a+i*h;
};

inline double Galerkin::phi(double x) // finite spline of the 1st order
{
    return (fabs(x)<1)? (1-fabs(x)) : 0;
};

inline double Galerkin::dPhi(double x) // derivative of phi(x)
{
    return (fabs(x)<1)? ((x<0)? 1 : -1) : 0;
};

inline double Galerkin::phi(long int i, double x) // finite spline for given mesh
{
    return phi((x-node(i))/h);
};

inline double Galerkin::dPhi(long int i, double x) // derivative of phi(i, x)
{
    return dPhi((x-node(i))/h)/h;
};

inline Galerkin::B::B(Galerkin& myP, int myC, int myR):
    c(myC),
    r(myR),
    p(myP)
{};

inline double Galerkin::B::operator()(double x)
{
    return p.eq.p(x)*p.dPhi(c,x)*p.dPhi(r,x)+
        p.eq.q(x)*p.phi(c,x)*p.phi(r,x);
};

inline Galerkin::FPhi::FPhi(Galerkin& myP, int myR):
    r(myR),
    p(myP)
{};

double Galerkin::FPhi::operator()(double x)
{
    return p.eq.f(x)*p.phi(r, x);
}

double Galerkin::distance(double* newu, double *oldu)
{
    double max;
    double temp;
    long int i;

    max = 0;
    for (i = 0; i < eq.n-1; i++)
        if ((temp = fabs(newu[i]-oldu[i])) > max)
            max = temp;
    return max;
}

void Galerkin::fillLinearSystem(
    LinearSystem* ls,
    long int n // n+1 = number of terms in the Galerkin's expansion 
)
{
    double x;
    double t1, t2, t3;
    long int i;

    h = (eq.b-eq.a)/n; // let's calculate the "smoothness" of splines

    if (fabs(eq.l0)>eq.eps) { // check if left edge condition is like
                  // y(a)=const
        t1 = integrate(eq.a, eq.a+h, B(*this, 0, 0))-
            eq.p(eq.a)*(eq.m0-1)/eq.l0;
        t2 = integrate(eq.a, eq.a+h, B(*this, 1, 0));
        t3 = integrate(eq.a, eq.a+h, FPhi(*this, 0));
        ls->k1 = -t2/t1;
        ls->n1 = t3/t1;
    } else {
        ls->k1 = 0;
        ls->n1 = eq.m0;
    }

    for (i = 0, x = eq.a; i < n-1; i++, x += h) {
        ls->a[i] = integrate(x, x+h, B(*this, i, i+1));
        ls->b[i] = integrate(x, x+2*h, B(*this, i+1, i+1));
        ls->c[i] = integrate(x+h, x+2*h, B(*this, i+2, i+1));
        ls->f[i] = integrate(x, x+2*h, FPhi(*this, i+1));
    }

    if (fabs(eq.l0)>eq.eps) { // check if right edge condition is like
                              // y(b)=const
        t1 = integrate(eq.b-h, eq.b, B(*this, n, n+1));
        t2 = integrate(eq.b-h, eq.b, B(*this, n+1, n+1))-
            eq.p(eq.b)*(eq.m1-1)/eq.l1;
        t3 = integrate(eq.b-h, eq.b, FPhi(*this, n+1));
        ls->k2 = -t1/t2;
        ls->n2 = t3/t2;
    } else {
        ls->k2 = 0;
        ls->n2 = eq.m1;
    }
}

//          n
//         ---
//         \    n    n
// y (x) = /   c  phi (x) (*) - Nth approximant of the solution
//  n      ---  j    j
//         j=0
//
double* Galerkin::solveDifEq(
    long int n // n+1 = number of terms in the Galerkin's expansion
)
{
    LinearSystem ls(n);
    double *y; // values of solution on the mesh
    double *c; // coefficients of the Galerkin's expansion
    long int i, j;
    double h; // granularity of the output mesh

    fillLinearSystem(&ls, n);
    c = ls.solve(); // now c has n+1 double elements
    y = new double[eq.n+1];
    h = (eq.b-eq.a)/eq.n;
    for (i = 0; i < eq.n+1; i++) { // for each node of the output mesh
        y[i] = 0;
        for (j = 0; j < n+1; j++) // calculation of (*)
            y[i] += c[j]*phi(j, eq.a+i*h);
    }
    delete [] c;
    return y;
}

void Galerkin::solve()
{
    long int m; // m+1 = number of terms in the Galerkin's expantion
    long int i;
    double *oldY;
    double *newY;
    double h; // granularity of the output mesh

    m = eq.n;
    oldY = NULL;
    newY = solveDifEq(m);
    do {
        if (oldY != NULL)
            delete [] oldY;
        oldY = newY;

// let's double the smoothness of the mesh
        m <<= 1;
        newY = solveDifEq(m);
    } while (distance(newY, oldY) >= eq.eps);

// Output of results
// It shouldn't be here!

    cout << "\nObtained results (number of terms in the Galerkin's expansion = "
        << m+1 << "):\n";
    h = (eq.b-eq.a)/eq.n;
    for (i = 0; i < eq.n+1; i++)
        cout << "y(" << eq.a+i*h << ") = " << newY[i] << "\n";
    cout << "Press any key...\n";
    char ch;
    cin.get(ch);

// End of output

    delete [] newY;
    delete [] oldY;
}