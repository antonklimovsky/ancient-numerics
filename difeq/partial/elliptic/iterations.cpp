//------------------------------------------------------------//
//  Method of interations for solving of the Neumann problem  //
//  for differential equation						          //
//                                                            //
//  Laplace(u)+px*u'+py*u'+q*u+f=0                            //
//                 x     y                                    //
//                                                            //
//	in the rectangle.                                         //
//                                                            //
//  Author: A. Klimovsky                                      //
//------------------------------------------------------------//
#include <iostream.h>
#include <math.h>
#include <conio.h>

#include "difeq.h"
#include "iterations.h"
#include "matrix.h"

Iterations::Iterations(DifEq& myEq):
	maxIterations(1000),
	scale(2),
	DifEqSolver(myEq)
{};

double Iterations::distance(Matrix& newu, Matrix& oldu)
{
	double max;
	double temp;
	long int i, j;

	max = 0;
	for (i = 0; i < n+1; i++)
		for (j = 0; j < n+1; j++)
			if ((temp = fabs(newu(i, j)-oldu(i, j))) > max)
				max = temp;
	return max;
}


//
// Calculation of the 1st mesh derivatives of y in the internal nodes
//

void Iterations::test(Matrix& a)
{
	long int i, j;
	char ch;

	clrscr();
	cout << "--------------\n";
	cout << "\nObtained results (number of nodes in mesh = "
		<< n+1 << ")\n";
	for (i = 0; i < n+1; i += gaps) {
		for (j = 0; j < n+1; j += gaps)
			printf("%5.3f ",  a(i, j));
		cout << endl;
	}
	cout << "--------------\n";
	cin.get(ch);
}

//
// Neumann boundary conditions
//
void Iterations::calculateBoundaryValues(Matrix& u)
{
	double t;
	long int i;

	for (i = 1, t = hX; i < n; i++, t += hX) { // bottom and top boundaries
		u(0, i) = (eq.bottom(t)*2*hY-4*u(1, i)+u(2, i))/(-3);
		u(n, i) = (eq.top(t)*2*hY+4*u(n-1, i)-u(n-2, i))/3;
	}
	for (i = 1, t = hY; i < n; i++, t += hY) { // left and right boundaries
		u(i, 0) = (eq.left(t)*2*hX-4*u(i, 1)+u(i, 2))/(-3);
		u(i, n) = (eq.right(t)*2*hX+4*u(i, n-1)-u(i, n-2))/3;
	}
}

/*
//
// Dirihlet boundary conditions
//
void Iterations::calculateBoundaryValues(Matrix& u)
{
	double t;
	long int i;

	for (i = 1, t = hX; i < n; i++, t += hX) { // bottom and top boundaries
		u(0, i) = eq.bottom(t);
		u(n, i) = eq.top(t);
	}
	for (i = 1, t = hY; i < n; i++, t += hY) { // left and right boundaries
		u(i, 0) = eq.left(t);
		u(i, n) = eq.right(t);
	}
}*/

void Iterations::meshDerivate(Matrix& newU, Matrix& oldU)
{
	double x, y;
	long int i, j;

	newU = oldU;
	for (i = 1, y = hY; i < n; i++, y += hY)
		for (j = 1, x = hX; j < n; j++, x += hX)
			newU(i, j) = oldU(i, j)+t*(
				(oldU(i-1, j)-2*oldU(i, j)+oldU(i+1, j))/(hY*hY)+
				(oldU(i, j-1)-2*oldU(i, j)+oldU(i, j+1))/(hX*hX)+
				eq.py(x, y)*(oldU(i+1, j)-oldU(i-1, j))/(2*hY)+
				eq.px(x, y)*(oldU(i, j+1)-oldU(i, j-1))/(2*hX)+
				eq.q(x, y)*oldU(i, j)+
				eq.f(x, y)
				);
}

//
// The iterations
//
void Iterations::solveDifEq()
{
	Matrix uPrev(n+1), uCur(n+1);
	long int i, j;
	long int counter; // of the iterations

	hX = eq.a/n; // let's calculate the "smoothness" of mesh
	hY = eq.b/n;

// let's calculate the iterations parameter
	t = 1/(2*sqr(sin(M_PI*hX/(2*eq.a)))/sqr(hX)+2*sqr(sin(M_PI*hY/(2*eq.b)))/sqr(hY)+
		   2*sqr(cos(M_PI*hX/(2*eq.a)))/sqr(hX)+2*sqr(cos(M_PI*hY/(2*eq.b)))/sqr(hY));
	cout << "t = " << t << endl;

	for (i = 0; i < n+1; i++) // the first approximant is 0
		for (j = 0; j < n+1; j++)
			uCur(i, j) = 0;

	counter = 0;
	do {
		uPrev = uCur;
		calculateBoundaryValues(uPrev);
//		test(uPrev);
		meshDerivate(uCur, uPrev); // let's make an iteration
//		cout << "a";
//		test(uCur);
		counter++;
		cout << "Iteration #" << counter << " residual:" <<
		fabs(distance(uCur, uPrev)) << "\n";
//	} while (distance(uCur, uPrev) > eq.eps && counter < maxIterations);
	} while (counter < maxIterations);

// Output of results
// It shouldn't be here!
	cout << "\nObtained results (number of nodes in mesh = "
		<< n+1 << ", " << counter << " iteration):\n";
	cout << "Residual = " << fabs(distance(uCur, uPrev)) << endl;
	for (i = 0; i < n+1; i += gaps) {
		for (j = 0; j < n+1; j += gaps)
			printf("%6.3f ",  uCur(i, j));
		cout << endl;
	}
	cout << "Press any key...\n";
	char ch;
	cin.get(ch);
// End of output
}

void Iterations::solve()
{
	long int i;

	n = eq.n * scale;
	gaps = 1 * scale;
	/*cout << "Please, input the iterations parameter: ";
	cin >> t;*/
	solveDifEq();
}

double Iterations::sqr(double x)
{
	return x*x;
}
