#include <iostream.h>
#include <math.h>

#include "difeq.h"
#include "iterations.h"

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
		eps = 1.0e-5;
		n = 20;
	};
	double p(double x) { return exp(sin(x)); };
	double g(double x, double y)
	{
		return sin(x)*exp(sin(x))*y+7*exp(x+sin(x));
	};
};

/*class MyDifEq: public DifEq {
public:
	MyDifEq()
	{
		a = 0;
		b = 1;
		m0 = 1;
		m1 = exp(1);
		l0 = 0;
		l1 = 0;
		eps = 1.0e-6;
		n = 20;
	};
	double p(double x) { return -exp(-x); };
	double g(double x, double y)
	{
		return 0;
	};
};*/

void main()
{
	DifEq* myDifEq;
	DifEqSolver* solver;

	myDifEq = new MyDifEq;
	solver = new Iterations(*myDifEq);
	solver->solve();
	delete solver;
	delete myDifEq;
}
