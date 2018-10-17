//-------------------------------------------------
// Author: Anton Klimovsky, root@ludus.kharkiv.com
//-------------------------------------------------
#ifndef __INTEGRATE_H
#define __INTEGRATE_H

//--------------------------------------------//
//                                  b         //
//                                  _         //
// Performs calculation of integral \ f(x)dx  //
//                                  -         //
//                                  a         //
// Calculation based on using of the mean     //
// point quadratures                          //
//                                            //
// This is really _very_ silly version!       //
//--------------------------------------------//
template<class Function> double myIntegrate(double a, double b, int n, Function f)
{
    double result;
    double h;
    double x;
    int i;

    h = (b-a)/n;
    result = 0;
    for (i = 0, x = a+h/2; i < n; i++, x+=h)
        result += f(x);
    return h*result;
}

template<class Function> double integrate(double a, double b, Function f)
{
    double result, oldResult;
    double eps = 1.0e-5;
    int n = 5;

    result = myIntegrate(a, b, n, f);
    do {
        oldResult = result;
        result = myIntegrate(a, b, n <<= 1, f);
    } while (fabs(oldResult-result)>eps);
    return result;
}

#endif
