#ifndef __GALERKIN_H
#define __GALERKIN_H

#include "difeq.h"
#include "matrix.h"

class Iterations: public DifEqSolver {
	const long int maxIterations;
	const long int scale;
	double hX;  // current granularity of the mesh along X axis
	double hY;  // current granularity of the mesh along Y axis
	long int n; // current nubmber of nodes in mesh
	double t;   // parameter of iterations
	long int gaps; // gaps between output nodes

	Matrix *pX, *pY, *q, *f, *tempU, *dX, *dY, *qU; // mesh functions

	void test(Matrix&);
	double distance(Matrix&, Matrix&);
	double sqr(double x);
	void calculateBoundaryValues(Matrix&);
	void solveDifEq();
	double applyStencil(int[3][3], Matrix&, long int, long int);
	double lambda1(Matrix&, long int, long int);
	double lambda2(Matrix&, long int, long int);
	double derive(double, double, double, double, double, double);
	double dx(Matrix&, long int, long int);
	double dy(Matrix&, long int, long int);
	double lambda12(Matrix&, long int, long int);
	double additionalFreeTerms(long int i, long int j);
	void precalculateMeshFunctions();
	void meshDerivate(Matrix&, Matrix&);
public:
	Iterations(DifEq&);
	void solve();
};

#endif
