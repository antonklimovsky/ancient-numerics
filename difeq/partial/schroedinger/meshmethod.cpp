//------------------------------------------------------------//
//  The mesh method for solving the Shroedinger equation
//
//  i*u = u"  + q*u' + f
//     t   xx      x
//  Author: A. Klimovsky
//------------------------------------------------------------//
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include <conio.h>

#include "difeq.h"
#include "meshmethod.h"
#include "progonka.h"

MeshMethod::MeshMethod(DifEq& myEq):
    DifEqSolver(myEq),
    scaleX(myEq.scaleX),
    scaleT(myEq.scaleT)
{};

inline double MeshMethod::sqr(double x)
{
    return x*x;
}

void MeshMethod::calculateNewLayer(cmplx* newU, cmplx* oldU, double t)
{
    using namespace Progonka;
    LinearSystem<cmplx> ls(n);
    long int i;
    cmplx im1(0, 1);
    double x;

    ls.k1 = 0;
    ls.n1 = eq.left(t);
/*  ls.k1 = 1;
    ls.n1 = -hT*eq.left(t);*/
    for (x = hX, i = 0; i < n-1; x += hX, i++) {
        ls.a[i] = ls.c[i] = -sigma*alpha;
        ls.b[i] = im1+2.0*sigma*alpha-hT*eq.q(x, t)/2;
        ls.f[i] = (1.0-sigma)*alpha*oldU[i] +
                  (im1-2.0*(1.0-sigma)*alpha+hT*eq.q(x, t-hT)/2)*oldU[i+1] +
                  (1.0-sigma)*alpha*oldU[i+2] +
                  hT*(eq.f(x, t)+eq.f(x, t-hT))/2.0;
    }
    ls.k1 = 0;
    ls.n1 = eq.right(t);
/*  ls.k1 = 1;
    ls.n1 = hT*eq.right(t);*/
    ls.solve();
    memcpy(newU, ls.x, sizeof(cmplx)*(n+1));
}

/*// Neumann problem
inline void MeshMethod::calculateBoundaryValues(cmplx* u, double t)
{
    u[0] = (eq.left(t)*2*hX-4*u[1]+u[2])/(-3);
    u[n] = (eq.right(t)*2*hX+4*u[n-1]-u[n-2])/3;
}*/

// Dirihlet problem
inline void MeshMethod::calculateBoundaryValues(cmplx* u, double t)
{
    u[0] = eq.left(t);
    u[n] = eq.right(t);
}

void MeshMethod::showLayer(cmplx* u)
{
    long int i;
    double sum;

    cout.setf(ios_base::fixed, ios_base::floatfield);
    cout.precision(3);

 /* for (i = 0; i < n+1; i += scaleX)
        cout <<  u[i];
    cout << endl;*/

    for (sum = i = 0; i < n+1; i++)
        sum += sqr(abs(u[i]))*hX;

    if (fabs(sum) > eq.eps)
        for (i = 0; i < n+1; i += scaleX)
            cout << sqr(abs(u[i]))*hX/sum << " ";

    cout << endl;
}

double MeshMethod::residual(cmplx* oldU, cmplx* newU, double t)
{
    long int i;
    cmplx im1(0, 1);
    double x;
    double result;
    double temp;

    result = 0;
    for (x = hX, i = 1; i < n; x += hX, i++) {
        temp = abs(im1*(newU[i]-oldU[i])/hT-
               (newU[i-1]-2.0*newU[i]+newU[i+1])/sqr(hX)-
               eq.q(x, t)*newU[i]-
               eq.f(x, t));
        if (temp > result)
            result = temp;
    }
//  cout << result << endl;
    return result;
}

//
// The mesh method
//
void MeshMethod::solveDifEq()
{
    long int i, j;
    double x, t;
    double eps;
    double temp;
    char ch;
    cmplx *oldU, *newU;
    cmplx im1(0, 1);

    oldU = new cmplx[n+1];
    newU = new cmplx[n+1];

    hX = eq.a/n; // let's calculate the granularity of the mesh
    hT = eq.T/m;

    alpha = hT/hX;
    sigma = 0.5;
    eps = 0;

    for (i = 0, x = 0; i < n+1; i++, x += hX)
        oldU[i] = eq.initial(x);
    calculateBoundaryValues(oldU, t);
    showLayer(oldU);

    for (i = 1, t = hT; i < m+1; i++, t += hT) {
        calculateNewLayer(newU, oldU, t);
        calculateBoundaryValues(newU, t);
        temp = residual(oldU, newU, t);
        if (temp > eps)
            eps = temp;
        memcpy(oldU, newU, sizeof(cmplx)*(n+1));
        if (!(i % scaleT))
            showLayer(oldU);
 //         cin.get(ch);
    }

    cout << "Residual " << eps << ". Press any key..." << endl;
    cin.get(ch);

    delete newU;
    delete oldU;
}

void MeshMethod::solve()
{
    long int i;

    n = eq.n * scaleX;
    m = eq.m * scaleT;
    solveDifEq();
}