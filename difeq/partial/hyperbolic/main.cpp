#include <iostream.h>
#include <math.h>

#include "difeq.h"
#include "meshmethod.h"

/*class MyDifEq: public DifEq {
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
	virtual double initialVelocity(double x)
	{
		return 0;
	}
	virtual double p(double x, double t)
	{
		return cos(x);
	};
	virtual double q(double x, double t)
	{
		return -sin(x);
	};
	virtual double f(double x, double t)
	{
		return -7*exp(x);
	};
};*/

class MyDifEq: public DifEq {
public:
	MyDifEq()
	{
		a = 1;
		T = 2;
		eps = 1.0e-3;
		n = 10;
		m = 20;
	};
	double left(double t) { return 0; };
	double right(double t) { return 0; };
	virtual double initial(double x)
	{
		return x*(1-x);
	}
	virtual double initialVelocity(double x)
	{
		return 0;
	}
	virtual double p(double x, double t)
	{
		return 0;
	};
	virtual double q(double x, double t)
	{
		return 0;
	};
	virtual double f(double x, double t)
	{
		return 0;
	};
};

/*class MyDifEq: public DifEq {
public:
	MyDifEq()
	{
		a = 1;
		T = 1;
		eps = 1.0e-3;
		n = 10;
		m = 10;
	};
	double left(double t) { return t; };
	double right(double t) { return a+t; };
	virtual double initial(double x)
	{
		return x;
	}
	virtual double initialVelocity(double x)
	{
		return 1;
	}
	virtual double p(double x, double t)
	{
		return 0;
	};
	virtual double q(double x, double t)
	{
		return 0;
	};
	virtual double f(double x, double t)
	{
		return 0;
	};
};*/

/*class MyDifEq: public DifEq {
public:
	MyDifEq()
	{
		a = 1;
		T = 1;
		eps = 1.0e-3;
		n = 10;
		m = 10;
	};
	double left(double t) { return 1; };
	double right(double t) { return 1; };
	virtual double initial(double x)
	{
		return x;
	}
	virtual double initialVelocity(double x)
	{
		return 1;
	}
	virtual double p(double x, double t)
	{
		return 0;
	};
	virtual double q(double x, double t)
	{
		return 0;
	};
	virtual double f(double x, double t)
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
