//---------------------------------------------------------//
// Calculation of the high-frequency Fourier coeffitients  //
// using the method of Filon.                              //
//---------------------------------------------------------//
#include <stdio.h>
#include <process.h>
#include <math.h>

//#define a 0.2
//#define b 1.2
double a = -M_PI;
double b = M_PI;

#define c 90
#define d 140
#define n0 5
#define eps 1.0e-6

double f(double x)
{
/*	if (x >= 0)
		return sqrt(x)*cos(x*x);
	else {
		puts("Error: imaginary result.");
		exit(1);
	}*/
	return sin(90*x);
}


double integrate(double omega, long int n)
{
	double h;
	double result;
	double temp, co, si, p, ReAp, ImAp;
	double x;
	double hdiv2;
	long int j;

	result = 0;
	hdiv2 = (h = (b-a)/n)/2;
	p = omega*h/2;
	if (fabs(p) < eps) {
		puts("Error: loss of accuracy is possible.");
		exit(1);
	}
	ReAp = sin(p)/p;
	ImAp = cos(p)-ReAp;

	for (j = 0, x = a; j < n-1; j++, x+=h) {
		temp = x+hdiv2;
		co = cos(omega*temp);
		si = sin(omega*temp);
		result += f(x)*(co*ImAp+si*ReAp)+
				  f(x+h)*(-co*ImAp+si*ReAp);
	}

	return h*result/2;
}

void main()
{
	double omega;
	double result, old_result;
	double t0 = (d-c)/n0;
	long int n;

	puts("Let's start the calculation...");
	for (omega = c; omega <= d; omega += t0) {
		n = 10*omega;
		result = integrate(omega, n);
		do {
			old_result = result;
			n <<= 1;
			result = integrate(omega, n);
		} while (fabs(result-old_result) >= eps);
		printf("omega = %lf, result = %lf\n",
					   omega, result);
	}
}
