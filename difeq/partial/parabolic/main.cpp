#include <iostream.h>
#include <math.h>

#include "difeq.h"
#include "meshmethod.h"

class MyDifEq: public DifEq {
public:
	MyDifEq()
	{
		a = 2;
		T = 2;
		eps = 1.0e-3;
		n = 10;
		m = 10;
	};
	double left(double t) { return 0; };
	double right(double t) { return 4; };
	virtual double initial(double x)
	{
		return 0;
	}
	virtual double p(double x)
	{
		return cos(x);
	};
	virtual double q(double x)
	{
		return -sin(x);
	};
	virtual double f(double x)
	{
		return -7*exp(x);
	};
};

/*class MyDifEq: public DifEq {
public:
	MyDifEq()
	{
		a = 1;
		T = 10;
		eps = 1.0e-3;
		n = 10;
		m = 10;
	};
	double left(double t) { return 1; };
	double right(double t) { return exp(1); };
	virtual double initial(double x)
	{
		return 0;
	}
	virtual double p(double x)
	{
		return 0;
	};
	virtual double q(double x)
	{
		return -1;
	};
	virtual double f(double x)
	{
		return 0;
	};
};*/

void main()
{
	DifEq* myDifEq;
	DifEqSolver* solver;

	myDifEq = new MyDifEq;
	solver = new MeshMethod(*myDifEq);
	solver->solve();
	delete solver;
	delete myDifEq;
}
