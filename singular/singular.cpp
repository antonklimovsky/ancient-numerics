//---------------------------------------------------//
// Calculates singular integrals with kernels of the //
// form 1/(x0-x).                                    //
//---------------------------------------------------//

#include <stdio.h>
#include <process.h>
#include <math.h>

#define a 0.2
#define b 1.2
#define n0 20
#define eps 1.0e-5
#define m 4

double c[] = {0.1739274225687284, 0.3260725774312716,
          0.3260725774312716, 0.1739274225687284};
double x[] = {-0.3611263115940442, -0.3399810435848646,
           0.3399810435848646,  0.3611263115940442};

// Integrable function
double f(double x)
{
    if (x >= 0)
        return sqrt(x)*cos(x*x);
    else {
        puts("Error: imaginary result.");
        exit(1);
    }
}

//------------------------------------------------------------------//
//                                b                                 //
//                                _                                 //
// Computes the singular integral \ f(x)dx/(x-x0)                   //
//                                -                                 //
//                                a                                 //
// Calculation uses the Gaussian quadratures of the                 //
// m-th order.                                                      //
// Parameter n specifies the granularity of the mesh used.          //
//------------------------------------------------------------------//
double integrate(double x0, int n)
{
    double h = (b-a)/n;
    double result;
    double temp;
    double z;
    int j, k;

    result = 0;
    for (j = 0; j < n; j++) {
        temp = 0;
        for (k = 0; k < m; k++) {
            z = a + (j+0.5)*h + h*x[k]/2;
            temp += c[k]*f(z)/(z-x0);
        }
        result += temp*h;
    }

    return result;
}

void main()
{
    double x0;
    double result, old_result;
    double h = (b-a)/n0;
    long int n;

    puts("Let's start the calculation...");
    for (x0 = a+h; x0 < b; x0 += h) {
        n = n0;
        result = integrate(x0, n);
        do {
            old_result = result;
            n <<= 1;
            result = integrate(x0, n);
        } while (fabs(result-old_result) >= eps);
        printf("x0 = %lf, result = %lf, after %i iterations\n", x0, result, n/n0);
    }
}
