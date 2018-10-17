#include <iostream.h>
#include <process.h>
#include "progonka.h"

LinearSystem::LinearSystem(int n0)
{
	n = n0;
	a = new double[n-1];
	b = new double[n-1];
	c = new double[n-1];
	f = new double[n-1];
}

LinearSystem::~LinearSystem()
{
	delete [] a;
	delete [] b;
	delete [] c;
	delete [] f;
}

double* LinearSystem::solveLinearSystem()
{
	double* alpha;
	double* beta;
	double* x;
	long int i;

	if(!(alpha = new double[n])) { cout << "Not enough memory"; exit(1); };
	if(!(beta = new double[n])) { cout << "Not enough memory"; exit(1); };
	x = new double[n+1];

	*alpha = k1;
	*beta = n1;
	for (i = 0; i < n-1; i++) {
		alpha[i+1] = -c[i]/(b[i]+a[i]*alpha[i]);
		beta[i+1] = (f[i]-a[i]*beta[i])/(b[i]+a[i]*alpha[i]);
	}

	x[n] = (k2*beta[n-1]+n2)/(1-k2*alpha[n-1]);
	for (i = n-1; i >= 0; i--)
		x[i] = alpha[i]*x[i+1]+beta[i];

	delete [] alpha;
	delete [] beta;
	return x;
}