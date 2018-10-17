//----------------------------------------------
// Solves folowing Fredholm's integral equation
// of the 2nd order:
//
// \begin{equation}
//     \label{eq:testFredholm2}
//     y(t)+\frac{1}{2}\int_0^1{e^{-ts}y(s)ds} = t+\frac{1-e^{-t}}{2}
// \end{equation}
//
// Author: A. Klimovsky, root@ludus.kharkiv.com
//----------------------------------------------
#include <iostream.h>
#include <math.h>

#include "lin_alg.h"
#include "gauss.h"
#include "integrator.h"

const double a = 0;
const double b = 1;
const double eps = 1.0e-3;

long int n0;
long int n = 50;

double h;
vector c;

double f(double t)
{
    return t+(1-exp(-t))/2;
}

double psi(int n, double t)
{
    int i;
    double temp;

    temp = 1;
    for (i = 1; i <= n; i++)
        temp *= t;
    return temp;
}

double phi(int n, double t)
{
    int i;
    double temp;

    temp = 0.5;
    for (i = 1; i <= n; i++)
        temp *= (-t)/(i);
    return temp;
}

class psiPhi {
    int i;
    int j;

public:
    psiPhi(int i0, int j0):
        i(i0),
        j(j0)
    {}

    double operator() (double t)
    {
        return psi(i, t)*phi(j, t);
    }
};

class psiF {
    int i;

public:
    psiF(int i0):
        i(i0)
    {}

    double operator() (double t)
    {
        return psi(i, t)*f(t);
    }
};

void inputMeshGranularity()
{
    cout << "Number of points along X axis: ";
    cin >> n0;
}

void solve()
{
    matrix A(n, n);
    vector f(n);
    long int i;
    long int j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            A[i][j] = integrate(a, b, psiPhi(i, j));
        A[i][i] += 1;
        f[i] = integrate(a, b, psiF(i));
    }
    c = calculate_inverse_matrix(A)*f;
}

double y(double t)
{
    int i;
    double temp;

    temp = 0;
    for (i = 0; i < n; i++)
        temp += c[i]*phi(i, t);

    return f(t)-temp;
}

void initialize()
{
    h = (b-a)/n0;
}

void output()
{
    long int i;
    double t;

    cout.setf(ios_base::fixed, ios_base::floatfield);
    cout.precision(3);

    for (i = 0, t = a; i < n0; i++, t += h)
        cout << "y(" << t << ") = " << y(t) << endl;
}

void main()
{
    inputMeshGranularity();
    initialize();
    solve();
    output();
}