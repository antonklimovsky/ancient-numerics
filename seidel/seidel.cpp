/*
	An implementation of the Seidel algorithm for solving the system of linear
	equations.
*/
#include <iostream.h>
#include <string.h>
#include <math.h>
#include <process.h>
#include "lin_alg.h"
//////////////////////////////////////////////////////////////////
// Algorithm starts here
//////////////////////////////////////////////////////////////////
#define EPS 1.0e-10
#define MAX_ITERATIONS 100
vector iterate(vector& x, matrix& a, vector& b)
{
	vector temp(b);
	int i, j;

	for (i = 0; i < b.get_size(); i++) {
		for (j = 0; j < i; j++)
			temp[i] -= a[i][j]*temp[j];

		for (j = i+1; j < b.get_size(); j++)
			temp[i] -= a[i][j]*x[j];

		if (fabs(a[i][i]) < EPS) {
			cout << "Can't process matrices with 0 on the main diagonal.\n";
			exit(1);
		}
		temp[i] *= 1/(a[i][i]);
	}
	return temp;
}
vector solve_linear_system(matrix& a, vector &b)
{
	int i;
	vector x_old(b), x_new(b.get_size());

	for (i = 0; i < MAX_ITERATIONS; i++) {
		x_new = iterate(x_old, a, b);
		if (abs(x_new-x_old) < EPS) {
			cout << i << " iterations have been done.\n";
			return x_new;
		}
		x_old = x_new;
	}
	cout << "No convergence after " << MAX_ITERATIONS << " steps.\n";
	exit(1);
}
//////////////////////////////////////////////////////////////////
// Main body of the program
//////////////////////////////////////////////////////////////////
void main()
{
	matrix a;
	vector b;
	vector x;
	int n;

	cout << "Input dimension of the system: ";
	cin >> n;

	cout << "Input matrix:\n";
	a = matrix(n, n);
	cin >> a;

	cout << "Input vector:\n";
	b = vector(n);
	cin >> b;

	x = solve_linear_system(a, b);
	cout << "Result:\n" << x << "\n";
	cout << "Product:\n" << a*x << "\n";
}
