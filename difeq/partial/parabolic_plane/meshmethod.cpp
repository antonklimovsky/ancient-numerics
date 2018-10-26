//------------------------------------------------------------//
//  The mesh method for solving of the parabolic equation     //
//                                                            //
//  u'=Laplace(u)+px*u'+py*u'+q*u+f                           //
//   t                x     y                                 //
//                                                            //
// in a rectangular domain with Neumann or Dirichle           //
// boundary conditions.                                       //
//                                                            //
// Author: A. Klimovsky                                       //
//------------------------------------------------------------//
#include <iostream.h>
#include <math.h>
#include <conio.h>

#include "difeq.h"
#include "meshmethod.h"
#include "matrix.h"

MeshMethod::MeshMethod(DifEq& myEq):
	scaleX(2),
	scaleT(1e+3),
	DifEqSolver(myEq)
{};

inline double MeshMethod::sqr(double x)
{
	return x*x;
}

//
// Second central derivative along X axis
//
inline double MeshMethod::lambda1(Matrix &u, long int i, long int j)
{
	return (u(i, j-1)-2*u(i, j)+u(i, j+1))/(hX*hX);
}

//
// Second central derivative along Y axis
//
inline double MeshMethod::lambda2(Matrix &u, long int i, long int j)
{
	return (u(i-1, j)-2*u(i, j)+u(i+1, j))/(hY*hY);
}

void MeshMethod::calculateNewLayer(Matrix& newU, Matrix& oldU)
{
	long int i, j;
	double x, y;

	for (i = 1, y = hY; i < n; i++, y += hY)
		for (j = 1, x = hX; j < n; j++, x+= hX)
			newU(i, j) = oldU(i, j)+hT*(
						lambda1(oldU, i, j)+lambda2(oldU, i, j)+
						eq.px(x, y)*(oldU(i+1, j)-oldU(i-1, j))/(2*hY) +
						eq.py(x, y)*(oldU(i, j+1)-oldU(i, j-1))/(2*hX) +
						eq.q(x, y)*oldU(i, j)+
						eq.f(x, y)
						);
}

// Neumann problem case
inline void MeshMethod::calculateBoundaryValues(Matrix& u, double t)
{
	long int i;
	double t1;

	for (i = 0, t1 = 0; i < n+1; i++, t1 += hY) { // bottom and top boundaries
		u(0, i) = (eq.bottom(t, t1)*2*hY-4*u(1, i)+u(2, i))/(-3);
		u(n, i) = (eq.top(t, t1)*2*hY+4*u(n-1, i)-u(n-2, i))/3;
	}
	for (i = 0, t1 = 0; i < n+1; i++, t1 += hX) { // left and right boundaries
		u(i, 0) = (eq.left(t, t1)*2*hX-4*u(i, 1)+u(i, 2))/(-3);
		u(i, n) = (eq.right(t, t1)*2*hX+4*u(i, n-1)-u(i, n-2))/3;
	}
}

/*// Dirihlet problem case
inline void MeshMethod::calculateBoundaryValues(Matrix& u, double t)
{
	double t1;
	long int i;

	for (i = 0, t1 = 0; i < n+1; i++, t1 += hX) { // bottom and top boundaries
		u(0, i) = eq.bottom(t, t1);
		u(n, i) = eq.top(t, t1);
	}
	for (i = 0, t1 = 0; i < n+1; i++, t1 += hY) { // left and right boundaries
		u(i, 0) = eq.left(t, t1);
		u(i, n) = eq.right(t, t1);
	}
}*/

void MeshMethod::showLayer(Matrix& u, double t)
{
 	long int i, j;
	char ch;

	clrscr();
	cout << "--------------\n";
	cout << "\nObtained results (t = " << t << ")\n";
	for (i = 0; i < n+1; i += scaleX) {
		for (j = 0; j < n+1; j += scaleX)
			printf("%6.3f ",  u(i, j));
		cout << endl;
	}

	cout << "--------------\n";
	cin.get(ch);
}

//
// The mesh method
//
void MeshMethod::solveDifEq()
{
	long int i, j;
	double x, t;
	char ch;
	Matrix oldU(n+1), newU(n+1);

	hX = eq.a/n; // let's calculate the "smoothness" of mesh
	hY = eq.b/k;
	hT = eq.T/m;

	cout << "Stability parameter (should be less then 0.5) = " << hT/sqr(hX) <<
		endl << "Press any key to start calculation..." << endl;
	cin.get(ch);

	for (i = 0; i < n+1; i++) // the first approximant is 0
		for (j = 0; j < n+1; j++)
			oldU(i, j) = 0;

	showLayer(oldU, 0);

	for (i = 1, t = hT; i < m+1; i++, t += hT) {
		calculateBoundaryValues(newU, t);
		calculateNewLayer(newU, oldU);
		oldU = newU;
		if (!(i % scaleT))
			showLayer(oldU, t);
	}

	cout << "That's all folks. Press any key..." << endl;
	cin.get(ch);
}

void MeshMethod::solve()
{
	long int i;

	n = eq.n * scaleX;
	k = eq.k * scaleX;
	m = eq.m * scaleT;
	solveDifEq();
}
