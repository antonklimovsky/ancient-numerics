#include <iostream.h>
#include <math.h>

#include "difeq.h"
#include "meshmethod.h"

double s[] = {1/10, 2/10, 3/10, 4/10};
double epsilon[] = {0.08, 0.06, 0.04};

class MyDifEq: public DifEq {

	inline double e(double x)
	{
		return (fabs(x)<=1) ? (1-fabs(x)) : 0;
	}

	inline double e(int k, int j, double x)
	{
		return e((x-a*s[k])/epsilon[j])/epsilon[j];
	}

	inline double p(double x)
	{
		return (fabs(x)<=1) ? 1 : 0;
	}

	inline double p(int k, int j, double x)
	{
		return p((x-a*s[k])/epsilon[j])/(2*epsilon[j]);
	}

public:
	double left(double t) { return 0; }

	double right(double t) { return 0; }

	double initial(double x)
	{
		return e(1, 1, x);
	}

	double q(double x, double t)
	{
		return p(2, 1, x);
	}

	double f(double x, double t)
	{
		return 0;
	}
	MyDifEq()
	{
		a = (2.5);
		T = (2.4);
		eps = (1.0e-3);
		n = (10);
		m = (10);
		scaleX = (25);
		scaleT = (12);
	}
};


void main()
{
	DifEq* eq;
	DifEqSolver* solver;

	eq = new MyDifEq;

	cout << "Number of points along X axis: ";
	cin >> (eq->n);
	cout << "Scale factor along X axis: ";
	cin >> (eq->scaleX);
	cout << "Number of points along T axis: ";
	cin >> (eq->m);
	cout << "Scale factor along T axis: ";
	cin >> (eq->scaleT);

	cout << "Ok, so X: " << (eq->scaleX*eq->n) << " Y: " << (eq->scaleT*eq->m)
		<< endl;

	solver = new MeshMethod(*eq);

	solver->solve();

	delete solver;
	delete eq;
	cin.get();
}