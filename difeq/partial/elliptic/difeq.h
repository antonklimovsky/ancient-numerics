#ifndef __DIFEQ_H
#define __DIFEQ_H

class DifEq {
public:

// Domain of definition [0;a]x[0;b] 

	double a;
	double b;

// Neumann edge conditions

	virtual double left(double y) = 0;	// x=0
	virtual double right(double y) = 0;	// x=a
	virtual double bottom(double x) = 0;// y=0
	virtual double top(double x) = 0;	// y=b

// Accuracy

	double eps;

// Number of nodes in output mesh

	long int n;

// Coefficients of equation -Laplace(u)=px*u'+py*u'+q*u-f
//                                          x     y

	virtual double px(double x, double y) = 0;
	virtual double py(double x, double y) = 0;
   	virtual double q(double x, double y) = 0;
   	virtual double f(double x, double y) = 0;
};

class DifEqSolver {
protected:
	DifEq& eq;
	long int n;
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