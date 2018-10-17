#ifndef __MESHMETHOD_H
#define __MESHMETHOD_H

#include "difeq.h"
#include <complex.h>

class MeshMethod: public DifEqSolver {
    typedef complex<double> cmplx;
    const long int scaleX; // scale factor between number of nodes in actual
                           // and output mesh
    const long int scaleT; // ----"----
    long int gaps; // gaps between output nodes
    double hX;     // current granularity of the mesh along X axis
    double hT;     // current granularity of the mesh along T axis
    long int n;    // nubmber of nodes in mesh along X axis
    long int m;    // nubmber of nodes in mesh along T axis
    cmplx sigma; // parameter
    cmplx alpha; // another parameter

    double sqr(double);
    void calculateNewLayer(cmplx*, cmplx*, double);
	void calculateBoundaryValues(cmplx*, double);
	double residual(cmplx*, cmplx*, double);
	void showLayer(cmplx*);
	void solveDifEq();
public:
    MeshMethod(DifEq&);
    void solve();
};

#endif
