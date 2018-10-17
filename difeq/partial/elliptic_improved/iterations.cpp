//------------------------------------------------------------//
//  The method of interations for solving of the Neumann problem
//  for differential equationof the form
//
//  Laplace(u)+px*u'+py*u'+q*u+f=0
//                 x     y
//
//  in a rectangular domain.
//
//  Author: A. Klimovsky
//------------------------------------------------------------//
#include <iostream.h>
#include <math.h>
#include <conio.h>

#include "difeq.h"
#include "iterations.h"
#include "matrix.h"

Iterations::Iterations(DifEq& myEq):
	maxIterations(10000),
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

	for (i = 0, t = 0; i < n+1; i++, t += hX) { // bottom and top boundaries
		u(0, i) = (eq.bottom(t)*hY-4*u(1, i)+3*u(2, i)-(4.0/3.0)*u(3, i)+u(4, i)/4.0)/(-25.0/12.0);
		u(n, i) = (eq.top(t)*hY+4*u(n-1, i)-3*u(n-2, i)+(4.0/3.0)*u(n-3, i)-u(n-4, i)/4.0)/(25.0/12.0);
	}
	for (i = 0, t = 0; i < n+1; i++, t += hY) { // left and right boundaries
		u(i, 0) = (eq.left(t)*hX-4*u(i, 1)+3*u(i, 2)-(4.0/3.0)*u(i, 3)+u(i, 4)/4.0)/(-25.0/12.0);
		u(i, n) = (eq.right(t)*hX+4*u(i, n-1)-3*u(i, n-2)+(4.0/3.0)*u(i, n-3)-u(i, n-4)/4.0)/(25.0/12.0);
	}
	u(0, 0) = (eq.bottom(0)*2*hY-4*u(1, 0)+u(2, 0))/(-3);
	u(n, 0) = (eq.top(0)*2*hY+4*u(n-1, 0)-u(n-2, 0))/3;
	u(0, n) = (eq.right(0)*2*hX+4*u(0, n-1)-u(0, n-2))/3;
	u(n, n) = (eq.right(hY)*2*hX+4*u(n, n-1)-u(n, n-2))/3;
}

/*//
// Dirihlet boundary conditions
//
void Iterations::calculateBoundaryValues(Matrix& u)
{
	double t;
	long int i;

	for (i = 0, t = 0; i < n+1; i++, t += hX) { // bottom and top boundaries
		u(0, i) = eq.bottom(t);
		u(n, i) = eq.top(t);
	}
	for (i = 0, t = 0; i < n+1; i++, t += hY) { // left and right boundaries
		u(i, 0) = eq.left(t);
		u(i, n) = eq.right(t);
	}
}*/

//
// Second central derivative along X axis
//
inline double Iterations::lambda1(Matrix &u, long int i, long int j)
{
	return (u(i, j-1)-2*u(i, j)+u(i, j+1))/(hX*hX);
}

//
// Second central derivative along Y axis
//
inline double Iterations::lambda2(Matrix &u, long int i, long int j)
{
	return (u(i-1, j)-2*u(i, j)+u(i+1, j))/(hY*hY);
}

/*inline double Iterations::dx(Matrix &u, long int i, long int j)
{
	return (u(i, j+1)-u(i, j-1))/(2*hX);
}

inline double Iterations::dy(Matrix &u, long int i, long int j)
{
	return (u(i+1, j)-u(i-1, j))/(2*hY);
}*/

inline double Iterations::derive
(
	double y0, double y1, double y2, double y3, double y4, double h
)
{
	return ((-25.0/12.0)*y0+4.0*y1-3.0*y2+(4.0/3.0)*y3-y4/4.0)/h;
}

inline double Iterations::dx(Matrix &u, long int i, long int j)
{
	if (j+4 < n+1)
		return derive(u(i, j), u(i, j+1), u(i, j+2), u(i, j+3), u(i, j+4), hX);
	else
		return derive(u(i, j), u(i, j-1), u(i, j-2), u(i, j-3), u(i, j-4), -hX);
}

inline double Iterations::dy(Matrix &u, long int i, long int j)
{
	if (i+4 < n+1)
		return derive(u(i, j), u(i+1, j), u(i+2, j), u(i+3, j), u(i+4, j), hX);
	else
		return derive(u(i, j), u(i-1, j), u(i-2, j), u(i-3, j), u(i-4, j), -hX);
}

double Iterations::applyStencil(int s[3][3], Matrix& u, long int myI, long int myJ)
{
	double result;
	long int i, j;

	result = 0;
	for (i = -1; i <= 1; i++)
		for (j = -1; j <= 1; j++)
			result += s[1+i][1+j]*u(myI+i, myJ+j);
	return  result;
}

//
// lambda12(u) = lambda1(lambda2(u))*hX*hY
//
inline double Iterations::lambda12(Matrix &u, long int myI, long int myJ)
{
	static int s[3][3] = { // stencil weights
		1, -2, 1,
		-2, 4, -2,
		1, -2, 1
	};
	return applyStencil(s, u, myI, myJ);
}

inline double Iterations::additionalFreeTerms(long int i, long int j)
{
	return (sqr(hX)*lambda1(*f, i, j)+sqr(hY)*lambda2(*f, i, j))/12;
}

void Iterations::meshDerivate(Matrix& newU, Matrix& oldU)
{
	long int i, j;

	newU = oldU;

	for (i = 0; i < n+1; i++)
		for (j = 1; j < n; j++) {
			(*dX)(i, j) = (*pX)(i, j)*dx(oldU, i, j);
			(*dY)(j, i) = (*pY)(j, i)*dy(oldU, j, i);
		}
	for (i = 0; i < n+1; i++) { // for problem of Neumann!
		(*dX)(i, 0) = (*pX)(i, 0)*eq.left(i*hY);
		(*dX)(i, n) = (*pX)(i, n)*eq.right(i*hY);
		(*dY)(0, i) = (*pY)(0, i)*eq.bottom(i*hX);
		(*dY)(n, i) = (*pY)(n, i)*eq.top(i*hX);
	}
	for (i = 0; i < n+1; i++)
		for (j = 0; j < n+1; j++)
			(*qU)(i, j) = oldU(i, j)*(*q)(i, j);

	for (i = 1; i < n; i++)
		for (j = 1; j < n; j++)
			newU(i, j) = oldU(i, j)+t*(
				lambda1(oldU, i, j) + lambda2(oldU, i, j) +
				(sqr(1/hX)+sqr(1/hY))*lambda12(oldU, i, j)/12 +
				sqr(hX)*(lambda1(*dX, i, j)+lambda1(*dY, i, j)+lambda1(*qU, i, j))/12+
				sqr(hY)*(lambda2(*dX, i, j)+lambda2(*dY, i, j)+lambda2(*qU, i, j))/12+
				(*dX)(i, j) +
				(*dY)(i, j) +
				(*qU)(i, j) +
				(*f)(i, j)+additionalFreeTerms(i, j)
				);
}

/*void Iterations::meshDerivate(Matrix& newU, Matrix& oldU)
{
	long int i, j;
	Matrix& intermediateU(*tempU);

	intermediateU = oldU;  // first step
	for (i = 1; i < n; i++)
		for (j = 1; j < n; j++)
			intermediateU(i, j) = oldU(i, j)+t*(
				lambda1(oldU, i, j)+
				(*pX)(i, j)*dx(oldU, i, j)+
				(*q)(i, j)*oldU(i, j)
			);

	newU = intermediateU;  // second step
	for (i = 1; i < n; i++)
		for (j = 1; j < n; j++)
			newU(i, j) = intermediateU(i, j)+t*(
				lambda2(intermediateU, i, j)+
				(sqr(1/hX)+sqr(1/hY))*lambda12(intermediateU, i, j)/12+
				(*pY)(i, j)*dy(intermediateU, i, j)+
				(*f)(i, j)+additionalFreeTerms(i, j)
			);
}*/

void Iterations::precalculateMeshFunctions()
{
	double x, y;
	long int i, j;

	for (i = 0, y = 0; i < n+1; i++, y += hY)
		for (j = 0, x = 0; j < n+1; j++, x += hX) {
			 (*pX)(i, j) = eq.px(x, y);
			 (*pY)(i, j) = eq.py(x, y);
			 (*q)(i, j) = eq.q(x, y);
			 (*f)(i, j) = eq.f(x, y);
		}
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
//	t = 5e-6;
	cout << "Iterations parameter t = " << t << endl;
	char ch;
	cin.get(ch);

	pX = new Matrix(n+1);
	pY = new Matrix(n+1);
	q = new Matrix(n+1);
	f = new Matrix(n+1);
	dX = new Matrix(n+1);
	dY = new Matrix(n+1);
	qU = new Matrix(n+1);
	tempU = new Matrix(n+1);
	precalculateMeshFunctions();

	for (i = 0; i < n+1; i++) // the first approximant is 0
		for (j = 0; j < n+1; j++)
			uCur(i, j) = 0;

	counter = 0;
	do {
		uPrev = uCur;
		calculateBoundaryValues(uPrev);
		meshDerivate(uCur, uPrev); // let's make an iteration
		counter++;
		cout << "Iteration #" << counter << " residual:" <<
			fabs(distance(uCur, uPrev)) << "\n";
//	} while (counter < maxIterations && fabs(distance(uCur, uPrev)) > eq.eps);
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
	cin.get(ch);

// End of output

	delete pX;
	delete pY;
	delete q;
	delete f;
	delete dX;
	delete dY;
	delete qU;
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

inline double Iterations::sqr(double x)
{
	return x*x;
}
