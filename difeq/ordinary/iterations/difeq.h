#ifndef __DIFEQ_H
#define __DIFEQ_H

class DifEq {
public:
// Domain of definition boundaries
	double a;
	double b;
// Edge conditions (y(a)+l0*y'(a)=m0) and (y(b)+l1*y'(b)=m1)
	double m0;
	double m1;
	double l0;
	double l1;
// Accuracy
	double eps;
// Number of nodes in output mesh
	long int n;
// Coefficients of non-linear equation -(p*y')'+g(x,y)=0
	virtual double p(double x) = 0;
	virtual double g(double x, double y) = 0;
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