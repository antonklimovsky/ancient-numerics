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
	void test(Matrix&);
	double distance(Matrix&, Matrix&);
	double sqr(double x);
	void calculateBoundaryValues(Matrix&);
	void solveDifEq();
	void meshDerivate(Matrix&, Matrix&);
public:
	Iterations(DifEq&);
	void solve();
};

#endif
