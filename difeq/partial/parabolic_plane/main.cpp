#include <iostream.h>
#include <math.h>

#include "difeq.h"
#include "meshmethod.h"

class MyDifEq: public DifEq {
public:
	MyDifEq()
	{
		a = 2;
		b = 2;
		T = 3;
		eps = 1.0e-3;
		n = 10;
		m = 10;
		k = 10;
	};
	double left(double t, double y) { return 0; };
	double right(double t, double y) { return 4; };
	double bottom(double t, double x) { return 0; };
	double top(double t, double x) { return 0; };
	virtual double px(double x, double y)
	{
		return cos(x);
	};
	virtual double py(double x, double y)
	{
		return 0;
	};
	virtual double q(double x, double y)
	{
		return -sin(x);
	};
	virtual double f(double x, double y)
	{
		return -7*exp(x);
	};
};

/*class MyDifEq: public DifEq {
public:
	MyDifEq()
	{
		a = 2;
		b = 2;
		T = 3;
		eps = 1.0e-3;
		n = 10;
		m = 10;
		k = 10;
	};
	double left(double t, double y) { return 1; };
	double right(double t, double y) { return 1; };
	double bottom(double t, double x) { return 0; };
	double top(double t, double x) { return 0; };
	virtual double px(double x, double y)
	{
		return 0;
	};
	virtual double py(double x, double y)
	{
		return 0;
	};
	virtual double q(double x, double y)
	{
		return 0;
	};
	virtual double f(double x, double y)
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
