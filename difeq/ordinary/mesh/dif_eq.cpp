//--------------------------------------
// Solves ODEs of the form y"+py'+qy=f
// with boundary conditions of the form
// y+l1*y'=m0), y+l2*y'=m1
//--------------------------------------

#include <iostream.h>
#include <math.h>
#include "progonka.h"



class DifEq {
public:
// Domain of definition boundaries
	double a;
	double b;
// Boundary conditions (y+l1*y'=m0) and (y+l2*y'=m1)
	double m0;
	double m1;
	double l1;
	double l2;
// Accuracy
	double eps;
// Number of nodes in output mesh
	long int n;
// Coefficients of equation y"+py'+qy=f
	virtual double p(double x) = 0;
	virtual double q(double x) = 0;
	virtual double f(double x) = 0;
};

class MyDifEq: public DifEq {
public:
	MyDifEq()
	{
		a = 0;
		b = 2;
		m0 = 0;
		m1 = 4.0e+300;
		l1 = -1.0e+300;
		l2 = 1.0e+300;
		eps = 1.0e-5;
		n = 20;
	};
	double p(double x) { return cos(x); };
	double q(double x) { return -sin(x); };
	double f(double x) { return 7*exp(x); };
};

/*class MyDifEq: public DifEq {
public:
	MyDifEq()
	{
		a = 0;
		b = 2;
		m0 = 0;
		m1 = 4;
		l1 = 0;
		l2 = 0;
		eps = 1.0e-5;
		n = 20;
	};
	double p(double x) { return cos(x); };
	double q(double x) { return -sin(x); };
	double f(double x) { return 7*exp(x); };
};*/

class TestDifEq: public DifEq {
public:
	TestDifEq()
	{
		a = 0;
		b = 1;
		m0 = 1;
		m1 = exp(1);
		l1 = 0;
		l2 = 0;
		eps = 1.0e-5;
		n = 20;
	};
	double p(double) { return 0; };
	double q(double) { return -1; };
	double f(double) { return 0; };
};

void discretizeEq(DifEq& eq, LinearSystem* ls, long int n)
{
	double h = (eq.b-eq.a)/n;
	double temp1, temp2;
	double x;
	long int i;

	temp2 = h*h;
	for (i = 0, x = eq.a+h; i < n-1; i++, x += h) {
		temp1 = eq.p(x)*h/2;
		ls->a[i] = (1-temp1)/temp2;
		ls->b[i] = -2/temp2+eq.q(x);
		ls->c[i] = (1+temp1)/temp2;
		ls->f[i] = eq.f(x);
	}
	ls->k1 = -(4*eq.l1*(ls->c[0])+(ls->b[0])*eq.l1)/
			 (temp1 = ((ls->a[0])*eq.l1+2*h*(ls->c[0])-3*(ls->c[0])*eq.l1));
	ls->n1 = (2*h*eq.m0*(ls->c[0])+(ls->f[0])*eq.l1)/temp1;

	ls->k2 = (4*eq.l2*(ls->a[n-2])+(ls->b[n-2])*eq.l2)/
			 (temp1 = ((-ls->c[n-2])*eq.l2+2*h*(ls->a[n-2])+3*(ls->a[n-2])*eq.l2));
	ls->n2 = (2*h*eq.m1*(ls->a[n-2])-(ls->f[n-2])*eq.l2)/temp1;
}

double distance(double* newu, double *oldu, long int newstep,
				long int newn)
{
	double max;
	double temp;
	long int i, j, oldstep;

	max = 0;
	oldstep = newstep >> 1;
	for (i = 0, j = 0; i < newn+1; i += newstep, j += oldstep)
		if ((temp = fabs(newu[i]-oldu[j])) > max)
			max = temp;
	return max;
}

double* solveDifEq(DifEq& eq, long int n)
{
	LinearSystem ls(n);

	discretizeEq(eq, &ls, n);
	return ls.solveLinearSystem();
}

void solveDifEqWithAccuracy(DifEq& eq)
{
	long int m;
	long int step;
	long int i;
	double *oldu;
	double *newu;

	m = eq.n;
	step = 1;
	oldu = NULL;
	newu = solveDifEq(eq, m);
	do {
		if (oldu != NULL)
			delete [] oldu;
		oldu = newu;
		m <<= 1;
		step <<= 1;
		newu = solveDifEq(eq, m);
	} while (distance(newu, oldu, step, m) >= eq.eps);
        
// Output of results
	cout << "\nObtained results (number of nodes in mesh = " << m << "):\n";
	for (i = 0; i < m+1; i += step)
		cout << "y(" << (eq.a+i*(eq.b-eq.a)/m) << ") = " << newu[i] << "\n";
	delete [] newu;
	delete [] oldu;
}

void main()
{
	MyDifEq myDifEq;
	TestDifEq testDifEq;
	char ch;

	solveDifEqWithAccuracy(myDifEq);
       	cout << "Press any key...\n";
	cin.get(ch);
	solveDifEqWithAccuracy(testDifEq);
	cout << "Press any key...\n";
	cin.get(ch);
}
