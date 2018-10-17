#ifndef __ITERATIONS_H
#define __ITERATIONS_H

#include "difeq.h"
#include "progonka.h"

class Iterations: public DifEqSolver {
	const long int maxIterations;
	const long int scale;
	double h;   // current granularity of the mesh
	long int n; // current nubmber of nodes in mesh
	double t;   // parameter of iterations
	long int gaps; // gaps between output nodes
	void test(double*, long int);
	double distance(double*, double *);
	double* derive(double*);
	double* derive2(double*);
	void fillLinearSystem(LinearSystem*, double*);
	double* solveDifEq();
public:
	void solve();
	Iterations(DifEq&);
};

#endif
