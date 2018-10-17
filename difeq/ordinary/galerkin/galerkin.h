#ifndef __GALERKIN_H
#define __GALERKIN_H

#include "difeq.h"
#include "progonka.h"
#include "integrator.h"

class Galerkin: public DifEqSolver {
	double h; // current granularity of the mesh
	double distance(double*, double *);
	void fillLinearSystem(LinearSystem*, long int);
	double* solveDifEq(long int);
	double node(long int i); // h should be already assigned!
	double phi(double x); // finite spline of the 1st order
	double dPhi(double x); // derivative of phi(x)
	double phi(long int i, double x); // finite spline for given mesh
	double dPhi(long int i, double x); // derivative of phi(i, x)
	class B { // for calculation of the linear system coefficients
		int c; // current col
		int r; // current row
		Galerkin& p; // reference to parent
	public:
		B(Galerkin& myP, int myC, int myR);
		double operator()(double x);
	};
	friend class B;
	class FPhi { // for calculation of the linear system free term
		int r; // current row
		Galerkin& p; // reference to parent
	public:
		FPhi(Galerkin& myP, int myR);
		double operator()(double x);
	};
	friend class FPhi;
public:
	Galerkin(DifEq& myEq);
	void solve();
};

#endif
