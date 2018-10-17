//------------------------------------------------------------//
//  The mesh method for solving of the hyperbolic equation    //
//                                                            //
//  u' =u" +p(x)u'+q(x)u+f                                    //
//   tt  xx      x                                            //
//                                                            //
//	in the rectangle.                                         //
//                                                            //
//  Author: A. Klimovsky                                      //
//------------------------------------------------------------//
#include <iostream.h>
#include <math.h>
#include <conio.h>

#include "difeq.h"
#include "meshmethod.h"

MeshMethod::MeshMethod(DifEq& myEq):
	scaleX(3),
	scaleT(1e+5),
	DifEqSolver(myEq)
{};

inline double MeshMethod::sqr(double x)
{
	return x*x;
}

void MeshMethod::calculateSecondLayer(double* oldU, double* veryOldU)
{
	long int i;
	double x;

	for (i = 1, x = hX; i < n; i++, x += hX)
		oldU[i] = veryOldU[i]+
			      hT*eq.initialVelocity(x)+
			      sqr(hT)*(
					(eq.initial(x+hX)-2*eq.initial(x+hX)+eq.initial(x-hX))/sqr(hX)+
					eq.p(x, 0)*(veryOldU[i+1]-veryOldU[i-1])/(2*hX)+
					eq.q(x, 0)*veryOldU[i]+
					eq.f(x, 0)
				  )/2;
}

void MeshMethod::calculateNewLayer
(
	double* newU, double* veryOldU, double* oldU, double t
)
{
	long int i;
	double x;

	for (i = 1, x = hX ; i < n; i++, x += hX)
		newU[i] = 2*oldU[i]-veryOldU[i]+sqr(hT)*(
				  (oldU[i+1]-2*oldU[i]+oldU[i-1])/sqr(hX)+
				  eq.p(x, t-hT)*(oldU[i+1]-oldU[i-1])/(2*hX)+
				  eq.q(x, t-hT)*oldU[i]+
                  eq.f(x, t-hT));
}

/*// Neumann problem boundary conditions
inline void MeshMethod::calculateBoundaryValues(double* u, double t)
{
	u[0] = (eq.left(t)*2*hX-4*u[1]+u[2])/(-3);
	u[n] = (eq.right(t)*2*hX+4*u[n-1]-u[n-2])/3;
}*/

// Dirihlet problem boundary conditions
inline void MeshMethod::calculateBoundaryValues(double* u, double t)
{
	u[0] = eq.left(t);
	u[n] = eq.right(t);
}

void MeshMethod::showLayer(double* u)
{
	long int i;

	for (i = 0; i < n+1; i += scaleX)
		printf("%6.3f ", u[i]);
	cout << endl;
}

//
// The mesh method
//
void MeshMethod::solveDifEq()
{
	long int i, j;
	double x, t;
	char ch;
	double *veryOldU, *oldU, *newU;

	veryOldU = new double[n+1];
	oldU = new double[n+1];
	newU = new double[n+1];

	hX = eq.a/n; // let's calculate the "smoothness" of mesh
	hT = eq.T/m;

	cout << "Stability parameter (should be less then 0.5) = " << hT/sqr(hX) <<
		endl << "Press any key to start calculation..." << endl;
	cin.get(ch);

	for (i = 0, x = 0; i < n+1; i++, x += hX)
		veryOldU[i] = eq.initial(x);

	showLayer(veryOldU);

	calculateSecondLayer(oldU, veryOldU);
	calculateBoundaryValues(oldU, hT);

	for (i = 2, t = 2*hT; i < m+1; i++, t += hT) {
		calculateNewLayer(newU, veryOldU, oldU, t);
        calculateBoundaryValues(newU, t);
		memcpy(veryOldU, oldU, sizeof(double)*(n+1));
		memcpy(oldU, newU, sizeof(double)*(n+1));
		if (!(i % scaleT))
			showLayer(oldU);
	}

	cout << "Press any key..." << endl;
	cin.get(ch);

	delete newU;
	delete oldU;
   	delete veryOldU;
}

void MeshMethod::solve()
{
	long int i;

	n = eq.n * scaleX;
	m = eq.m * scaleT;
	solveDifEq();
}