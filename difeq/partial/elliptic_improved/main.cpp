#include <iostream.h>
#include <math.h>

#include "difeq.h"
#include "iterations.h"

class MyDifEq: public DifEq {
public:
	MyDifEq()
	{
		a = 2;
		b = 2;
		eps = 1.0e-3;
		n = 10;
	};
	double left(double y) { return 0; };
	double right(double y) { return 4; };
	double bottom(double x) { return 0; };
	double top(double x) { return 0; };
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
		a = 1;
		b = 1;
		eps = 1.0e-3;
		n = 10;
	};
	double left(double y) { return 1; };
	double right(double y) { return exp(1); };
	double bottom(double x) { return 0; };
	double top(double x) { return 0; };
	virtual double px(double x, double y)
	{
		return 1;
	};
	virtual double py(double x, double y)
	{
		return 0;
	};
	virtual double q(double x, double y)
	{
		return -1;
	};
	virtual double f(double x, double y)
	{
		return -exp(x);
	};
};*/

/*class MyDifEq: public DifEq { // solution = x^2+C
public:
	MyDifEq()
	{
		a = 1;
		b = 1;
		eps = 1.0e-5;
		n = 10;
	};
	double left(double y) { return 0; };
	double right(double y) { return 2; };
	double bottom(double x) { return 0; };
	double top(double x) { return 0; };
	virtual double px(double x, double y)
	{
		return 1;
	};
	virtual double py(double x, double y)
	{
		return 0;
	};
	virtual double q(double x, double y)
	{
		return 1;
	};
	virtual double f(double x, double y)
	{
		return -(2+2*x+x*x);
	};
};*/

/*class MyDifEq: public DifEq { // solution = x^2+C
public:
	MyDifEq()
	{
		a = 1;
		b = 1;
		eps = 1.0e-4;
		n = 10;
	};
	double left(double y) { return 0; };
	double right(double y) { return 2; };
	double bottom(double x) { return 0; };
	double top(double x) { return 0; };
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
		return -2;
	};
};*/

/*class MyDifEq: public DifEq { // solution = x^2
public:
	MyDifEq()
	{
		a = 1;
		b = 1;
		eps = 1.0e-4;
		n = 10;
	};
	double left(double y) { return 0; };
	double right(double y) { return 1; };
	double bottom(double x) { return x*x; };
	double top(double x) { return x*x; };
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
		return -2;
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
