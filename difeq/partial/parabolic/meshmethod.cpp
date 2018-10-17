//------------------------------------------------------------//
//  The mesh method for solving of the parabolic equation     //
//                                                            //
//  u'=u" +p(x)u'+q(x)u+f                                     //
//   t  xx      x                                             //
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

void MeshMethod::calculateNewLayer(double* newU, double* oldU)
{
	long int i;
	double x;

	for (i = 1, x = hX ; i < n; i++, x += hX)
		newU[i] = oldU[i]+hT*(
		          (oldU[i+1]-2*oldU[i]+oldU[i-1])/sqr(hX)+
				  eq.p(x)*(oldU[i+1]-oldU[i-1])/(2*hX)+
				  eq.q(x)*oldU[i]+
                  eq.f(x));
}

// Neumann problem case
inline void MeshMethod::calculateBoundaryValues(double* u, double t)
{
	u[0] = (eq.left(t)*2*hX-4*u[1]+u[2])/(-3);
	u[n] = (eq.right(t)*2*hX+4*u[n-1]-u[n-2])/3;
}

/*// Dirihlet problem case
inline void MeshMethod::calculateBoundaryValues(double* u, double t)
{
	u[0] = eq.left(t);
	u[n] = eq.right(t);
}*/

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
	double *oldU, *newU;

	oldU = new double[n+1];
	newU = new double[n+1];

	hX = eq.a/n; // let's calculate the "smoothness" of mesh
	hT = eq.T/m;

	cout << "Stability parameter (should be less then 0.5) = " << hT/sqr(hX) <<
		endl << "Press any key to start calculation..." << endl;
	cin.get(ch);

	for (i = 0, x = 0; i < n+1; i++, x += hX)
		oldU[i] = eq.initial(x);
	calculateBoundaryValues(oldU, t);		
	showLayer(oldU);

	for (i = 1, t = hT; i < m+1; i++, t += hT) {
		calculateBoundaryValues(newU, t);
		calculateNewLayer(newU, oldU);
		memcpy(oldU, newU, sizeof(double)*(n+1));
		if (!(i % scaleT))
			showLayer(oldU);
	}

	cout << "Press any key..." << endl;
	cin.get(ch);

	delete newU;
	delete oldU;
}

void MeshMethod::solve()
{
	long int i;

	n = eq.n * scaleX;
	m = eq.m * scaleT;
	solveDifEq();
}