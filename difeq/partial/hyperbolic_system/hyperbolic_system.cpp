#include <iostream.h>
#include <math.h>
#include <mem.h>

const double l = 1.3;
const double a[] = {2, -2};
const double hT0 = 0.35;

double* oldU1;
double* oldU2;
double* newU1;
double* newU2;
double hX, hT;
double t;
double gamma;
double tMax;
long int n;
long int m;
long int n0;
long int m0;
long int scaleX;
long int scaleT;

inline int theta(double x) // theta(-x) = Heviside's function
{
    return (x < 0);   
}

// potential matrix elements

double q11(double x, double t)
{
    return log(1+fabs(x));
}

double q21(double x, double t)
{
    return cos(x);
}

double q12(double x, double t)
{
    return sin(x);
}

double q22(double x, double t)
{
    return 0;
}

// free term

double f1(double x, double t) 
{
    return t;
}

double f2(double x, double t)
{
    return 0;
}

// initial condition (t = 0)

double phi1(double x)
{
    return cos(x);
}

double phi2(double x)
{
    return x*exp(-x);
}

// boundary condition

double psi1(double t)
{
    return phi1(l*theta(a[0]));
}

double psi2(double t)
{
    return phi2(l*theta(a[1]));
}


void initialLayer()
{
    long int i;
    double x;

    for (i = 0, x = 0; i < n; i++, x += hX) {
        oldU1[i] = phi1(x);
        oldU2[i] = phi2(x);
    }
}

void calcNextLayer()
{
    long int i;
    double x;

    if (a[0] > 0) {
        newU1[0] = psi1(t);
        for (i = 1, x = hX; i < n; i++, x += hX)
            newU1[i] = (a[0]*gamma*newU1[i-1]+
                       hT*(q11(x, t)*oldU1[i]+q12(x, t)*oldU2[i])+
                       oldU1[i]+
                       hT*f1(x, t))/(1+a[0]*gamma);
    } else {
        newU1[n-1] = psi1(t);
        for (i = n-2, x = l-hX; i >= 0; i--, x -= hX)
            newU1[i] = (-a[0]*gamma*newU1[i+1]+
                       hT*(q11(x, t)*oldU1[i]+q12(x, t)*oldU2[i])+
                       oldU1[i]+
                       hT*f1(x, t))/(1-a[0]*gamma);
    }

    if (a[1] > 0) {
        newU2[0] = psi2(t);
        for (i = 1, x = hX; i < n; i++, x += hX)
            newU2[i] = (a[1]*gamma*newU2[i-1]+
                       hT*(q21(x, t)*oldU2[i]+q22(x, t)*oldU2[i])+
                       oldU2[i]+
                       hT*f2(x, t))/(1+a[1]*gamma);
    } else {
        newU2[n-1] = psi2(t);
        for (i = n-2, x = l-hX; i >= 0; i--, x -= hX)
            newU2[i] = (-a[1]*gamma*newU2[i+1]+
                       hT*(q21(x, t)*oldU2[i]+q22(x, t)*oldU2[i])+
                       oldU2[i]+
                       hT*f2(x, t))/(1-a[1]*gamma);
    }
}

void inputMeshGranularity()
{
    cout << "Number of points along X axis: ";
    cin >> n0;
    cout << "Scale factor along X axis: ";
    cin >> scaleX;
    n = n0*scaleX;

    cout << "Number of points along T axis: ";
    cin >> m0;
    tMax = hT0*m0;
    cout << "Scale factor along T axis: ";
    cin >> scaleT;
    m = m0*scaleT;

    cout << "Ok, so X: " << n << " T: " << m << endl;
}

void initialize()
{
    oldU1 = new double[n];
    oldU2 = new double[n];
    newU1 = new double[n];
    newU2 = new double[n];

    hX = l/n;
    hT = tMax/m;
    gamma = hT/hX;
}

void destroy()
{
    delete [] oldU1;
    delete [] oldU2;
    delete [] newU1;
    delete [] newU2;
}

void oldify()
{
    memcpy(oldU1, newU1, sizeof(double)*n);
    memcpy(oldU2, newU2, sizeof(double)*n);
}

void outputLayer()
{
    long int i;

    cout.setf(ios_base::fixed, ios_base::floatfield);
    cout.precision(3);

    for (i = 0; i < n; i += scaleX)
        cout << oldU1[i] << " ";
    cout << endl;
    for (i = 0; i < n; i += scaleX)
        cout << oldU2[i] << " ";
    cout << endl << endl;
}

void calcResidual()
{
    double residual;
    double temp;
    double x;
    long int i;

    residual = 0;
    for (i = 0, x = 0; i < n-1; i++, x += hX) {
        temp = fabs((newU1[i] - oldU1[i])/hT+
                     a[0]*(oldU1[i+1]-oldU1[i])/hX-
                     q11(x, t)*oldU1[i]+q12(x, t)*oldU2[i]-
                     f1(x, t)
                     );
        temp += fabs((newU2[i] - oldU2[i])/hT+
                     a[1]*(oldU2[i+1]-oldU2[i])/hX-
                     q21(x, t)*oldU1[i]+q22(x, t)*oldU2[i]-
                     f2(x, t)
                     );
        if (residual < temp)
            residual = temp;
    }
    cout << "Residual on the last layer: " << residual << endl;
}

void solve()
{
    long int i;

    inputMeshGranularity();
    initialize();
    initialLayer();
    outputLayer();
    for (i = 1, t = 0; i < m; i++, t += hT) {
        calcNextLayer();
        if (i == m-1)
            calcResidual();
        oldify();
        if (!(i % scaleT))
            outputLayer();
    }
    destroy();
}

void main()
{
    char ch;

    solve();
    
    cin.get(ch);
}