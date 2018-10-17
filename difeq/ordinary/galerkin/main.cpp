#include <iostream.h>
#include <math.h>

#include "difeq.h"
#include "galerkin.h"
#include "integrator.h"

class MyDifEq: public DifEq {
public:
	MyDifEq()
	{
		a = 0;
		b = 2;
		m0 = 0;
		m1 = 4;
		l0 = 0;
		l1 = 0;
		eps = 1.0e-3;
		n =20;
	};
	double p(double x) { return exp(sin(x)); };
	double q(double x) { return sin(x)*exp(sin(x)); };
	double f(double x) { return -7*exp(x+sin(x)); };
};

/*class MyDifEq: public DifEq {
public:
	MyDifEq()
	{
		a = 0;
		b = 2;
		m0 = 0;
		m1 = 4.0e+300;
		l0 = -1.0e+300;
		l1 = 1.0e+300;
		eps = 1.0e-3;
		n =20;
	};
	double p(double x) { return exp(sin(x)); };
	double q(double x) { return sin(x)*exp(sin(x)); };
	double f(double x) { return -7*exp(x+sin(x)); };
};*/

void main()
{
	DifEq* myDifEq;
	DifEqSolver* solver;

	myDifEq = new MyDifEq;
	solver = new Galerkin(*myDifEq);
	solver->solve();
	delete solver;
	delete myDifEq;
}
