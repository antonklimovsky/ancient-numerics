#ifndef __DIFEQ_H
#define __DIFEQ_H

class DifEq {
public:

// Domain of definition [0;a]x[0;T] 

	double a;
	double T;

// Neumann edge conditions

	virtual double left(double t) = 0;	 // x=0
	virtual double right(double t) = 0;	 // x=a
	virtual double initial(double x) = 0;// t=0

// Accuracy

	double eps;

// Number of nodes in output mesh

	long int n, m;

// Coefficients of equation -Laplace(u)=px*u'+py*u'+q*u-f
//                                          x     y

	virtual double p(double x) = 0;
	virtual double q(double x) = 0;
   	virtual double f(double x) = 0;
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