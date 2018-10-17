//--------------------------------------------------------------//
//  The Method of interations for solving a boundary problem for//
//  the 2nd order ODEs.                                         //
//                                                              //
//  Author: A. Klimovsky                                        //
//--------------------------------------------------------------//
#include <iostream.h>
#include <math.h>

#include "difeq.h"
#include "progonka.h"
#include "iterations.h"

Iterations::Iterations(DifEq& myEq):
	maxIterations(500),
	scale(5),
	DifEqSolver(myEq)
{};

double Iterations::distance(double* newu, double *oldu)
{
	double max;
	double temp;
	long int i;

	max = 0;
	for (i = 0; i < n; i++)
		if ((temp = fabs(newu[i]-oldu[i])) > max)
			max = temp;
	return max;
}

//
// Calculation of the 1st mesh derivatives of y in the internal nodes
//
double* Iterations::derive(double* y)
{
	double* dy = new double[n-1];
	long int i;

	for(i = 1; i < n; i++)
		dy[i-1] = (y[i+1]-y[i-1])/(2*h);

	return dy;
}

//
// Calculation of the 2nd mesh derivatives of y in the internal nodes
//
double* Iterations::derive2(double* y)
{
	double* dy = new double[n-1];
	long int i;

	for(i = 1; i < n; i++)
		dy[i-1] = (y[i-1]-2*y[i]+y[i+1])/(h*h);

	return dy;
}

void Iterations::test(double *a, long int n)
{
	long int i;
	char ch;

	cout << "--------------\n";
	for (i = 0; i < n; i++)
		cout << a[i] << endl;
	cout << "--------------\n";
	cin.get(ch);
}

void Iterations::fillLinearSystem(
	LinearSystem* ls,
	double* yPrev
)
{
	double temp1, temp2;
	double x;
	long int i;
	double* yDerivation; // temporary vars for calc. of the linear system free term
	double* yDerivation2;
	double* p;
	double* pDerivation;

	// Calculation of the free term
	yDerivation = derive(yPrev);
	yDerivation2 = derive2(yPrev);
	p = new double[n+1];
	for (i = 0, x = eq.a; i < n+1; i++, x += h)
		p[i] = eq.p(x);
	pDerivation = derive(p);
	for (i = 0, x = eq.a+h; i < n-1; i++, x += h)
		ls->f[i] = yDerivation2[i]+t*(-pDerivation[i]*yDerivation[i]-
			eq.p(x)*yDerivation2[i]+eq.g(x, yPrev[i+1]));
	delete [] yDerivation;
	delete [] yDerivation2;
	delete [] pDerivation;
	delete [] p;

	temp2 = 1/(h*h);
	for (i = 0, x = eq.a+h; i < n-1; i++, x += h) {
		ls->a[i] = temp2;
		ls->b[i] = -2*temp2;
		ls->c[i] = temp2;
	}

	ls->k1 = -(4*eq.l0*(ls->c[0])+(ls->b[0])*eq.l0)/
			 (temp1 = ((ls->a[0])*eq.l0+2*h*(ls->c[0])-3*(ls->c[0])*eq.l0));
	ls->n1 = (2*h*eq.m0*(ls->c[0])+(ls->f[0])*eq.l0)/temp1;

	ls->k2 = (4*eq.l1*(ls->a[n-2])+(ls->b[n-2])*eq.l1)/
			 (temp1 = ((-ls->c[n-2])*eq.l1+2*h*(ls->a[n-2])+3*(ls->a[n-2])*eq.l1));
	ls->n2 = (2*h*eq.m1*(ls->a[n-2])-(ls->f[n-2])*eq.l1)/temp1;
}

//
// The iterations
//
double* Iterations::solveDifEq()
{
	LinearSystem ls(n);
	double *yPrev, *yCur; 
	long int i;
	long int counter; // of the iterations

	h = (eq.b-eq.a)/n; // let's calculate the "smoothness" of mesh
	yCur = new double[n+1];
	yPrev = new double[1];
	for (i = 0; i < n+1; i++) // the first approximant is 0
		yCur[i] = 0;

	counter = 0;
	do {
		delete yPrev;
		yPrev = yCur;
		fillLinearSystem(&ls, yPrev); // let's prepare inversion
		yCur = ls.solve(); // let's invert our mesh differential operator
		counter++;
	} while (distance(yCur, yPrev) > eq.eps && counter < maxIterations);

// Output of results
// FIXME: It shouldn't be here!
	cout << "\nObtained results (number of nodes in mesh = "
		<< n+1 << ", " << counter << " iteration):\n";
	for (i = 0; i < n+1; i += gaps)
		cout << "y(" << eq.a+i*h << ") = " << yCur[i] << "\n";
	cout << "Press any key...\n";
	char ch;
	cin.get(ch);
// End of output

	delete yPrev;
	return yCur;
}

void Iterations::solve()
{
	long int i;
	double *y;

	n = eq.n << scale;
	gaps = 1 << scale;
	/*cout << "Please, input the iterations parameter: ";
	cin >> t;*/
	t = 0.05;

	y = solveDifEq();

	delete [] y;
}
