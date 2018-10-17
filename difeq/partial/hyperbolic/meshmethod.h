#ifndef __GALERKIN_H
#define __GALERKIN_H

#include "difeq.h"

class MeshMethod: public DifEqSolver {
	const long int scaleX;// scale factor between number of nodes in actual
						  // and output mesh
    const long int scaleT;
   	long int gaps; // gaps between output nodes
	double hX;  // current granularity of the mesh along X axis
	double hT;  // current granularity of the mesh along T axis
	long int n; // nubmber of nodes in mesh along X axis
	long int m; // nubmber of nodes in mesh along T axis

	double sqr(double);
	void calculateSecondLayer(double*, double*);
	void calculateNewLayer(double* newU, double* veryOldU, double* oldU, double t);
	void calculateBoundaryValues(double*, double);
	void showLayer(double*);
	void solveDifEq();
public:
	MeshMethod(DifEq&);
	void solve();
};

#endif
