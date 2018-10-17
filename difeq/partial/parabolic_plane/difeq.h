#ifndef __DIFEQ_H
#define __DIFEQ_H

class DifEq {
public:

// Domain of definition [0;a]x[0;T] 

	double a;
	double b;
	double T;

// Edge conditions

	virtual double left(double t, double y) = 0;	 // x=0
	virtual double right(double t, double y) = 0;	 // x=a
	virtual double bottom(double t, double x) = 0;
	virtual double top(double t, double x) = 0;

// Accuracy

	double eps;

// Number of nodes in output mesh

	long int n, m, k;

// Coefficients of equation u'=Laplace(u)+px*u'+py*u'+q*u+f
//                           t                x     y

	virtual double px(double x, double y) = 0;
	virtual double py(double x, double y) = 0;
	virtual double q(double x, double y) = 0;
	virtual double f(double x, double y) = 0;
};

class DifEqSolver {
protected:
	DifEq& eq;
public:
        DifEqSolver(DifEq& myEq):
                eq(myEq)
        {};
        void setEquation(DifEq& myEq)
        {
                eq = myEq;
        };
        virtual void solve() = 0;
};

#endif