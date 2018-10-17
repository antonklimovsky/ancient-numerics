#ifndef __GALERKIN_H
#define __GALERKIN_H

#include "difeq.h"
#include "matrix.h"

class MeshMethod: public DifEqSolver {
	const long int scaleX;// scale factor between number of nodes in actual
						  // and output mesh
    const long int scaleT;
   	long int gaps; // gaps between output nodes
	double hX;  // current granularity of the mesh along X axis
	double hY;  // current granularity of the mesh along Y axis
	double hT;  // current granularity of the mesh along T axis
	long int n; // nubmber of nodes in mesh along X axis
	long int m; // nubmber of nodes in mesh along T axis
	long int k;

	double sqr(double);
	double lambda1(Matrix&, long int, long int);
	double lambda2(Matrix&, long int, long int);
	void calculateNewLayer(Matrix&, Matrix&);
	void calculateBoundaryValues(Matrix&, double);
	void showLayer(Matrix&, double);
	void solveDifEq();
public:
	MeshMethod(DifEq&);
	void solve();
};

#endif
