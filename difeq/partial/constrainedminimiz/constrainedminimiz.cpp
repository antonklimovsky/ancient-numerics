//-------------------------------------------------
// Multidimesional constrained minimization.
//
// Steepest descent method with outer penalties.
//
// Author: Anton Klimovsky, root@ludus.kharkiv.com
//-------------------------------------------------
#include <iostream.h>
#include <math.h>

const double eps = 1e-4;
const double tBoundary = 20;
const double T = 1.0e+6;
const double maxIterations = 1000;

inline double sqr(double x)
{
    return x*x;
}

double h1(double x, double y)
{
    return sqr(sqr(x-2))+sqr(sqr(y-3))-16;
}

double h2(double x, double y)
{
    return sqr(x-2)+2-y;
}

double h3(double x, double y)
{
    return y-4;
}

double f(double x, double y)
{
    return sqr(x)+5*sqr(y);
}

inline double truncate(double x)
{
    return (x > 0)? sqr(x) : 0;
}

double penalty(double x, double y, double T)
{
    return T*(truncate(h1(x, y))+
              truncate(h2(x, y))+
              truncate(h3(x, y)));
}

double f0(double x, double y)
{
    return f(x, y)+penalty(x, y, T);
}

class G {
    double x;
    double y;
    double dirX;
    double dirY;

public:
    G(double x0, double y0, double dirX0, double dirY0):
        x(x0),
        y(y0),
        dirX(dirX0),
        dirY(dirY0)
    {}
    double operator() (double t)
    {
        return f0(x-dirX*t, y-dirY*t);
    }
};

template<class F> void dihotomy(
    double a,
    double b,
    F f,
    double eps,
    double* x0
)
{
    double c;

    while ((b-a) > eps) {
        c = (a+b)/2;
        if (f((a+c)/2) < f((c+b)/2))
            b = c;
        else
            a = c;
    }
    *x0 = c;
}

void calc()
{
    double oldX;
    double oldY;
    double newX;
    double newY;
    double t0;
    double gradX;
    double gradY;
    double norm;
    long int i;

    newX = 2;
    newY = 3;
    i = 0;
    do {
        i++;

        oldX = newX;
        oldY = newY;

        gradX = (f0(oldX+eps, oldY)-f0(oldX-eps, oldY))/(2*eps);
        gradY = (f0(oldX, oldY+eps)-f0(oldX, oldY-eps))/(2*eps);

        norm = sqrt(sqr(gradX)+sqr(gradY));

        gradX /= norm;
        gradY /= norm;

        dihotomy(0, tBoundary, G(oldX, oldY, gradX, gradY), eps, &t0);

        newX = oldX-t0*gradX;
        newY = oldY-t0*gradY;
    }
    while (fabs(newX-oldX)+fabs(newY-oldY) > eps && i < maxIterations);

    cout << i << " iterations done..." << endl;

    cout << "x: " << newX << ", y: " << newY <<
         ", f(x, y): " << f(newX, newY) <<
         ", penalty(x, y): " << penalty(newX, newY, T) << endl;

    gradX = (f0(newX+eps, oldY)-f0(newX-eps, newY))/(2*eps);
    gradY = (f0(newX, newY+eps)-f0(newX, newY-eps))/(2*eps);
    cout << "gradX: " << gradX << " gradY: " << gradY << endl;
}

void main()
{
    calc();    
}