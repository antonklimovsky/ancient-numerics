/*
	An implementation of the Square Roots algorithm for solving the systems of
	linear equations with symmetric matrices.
*/
#include <iostream.h>
#include <string.h>
#include <math.h>
#include <process.h>
#include "lin_alg.h"
//////////////////////////////////////////////////////////////////
// Algorithm starts here
//////////////////////////////////////////////////////////////////
vector solve_linear_system(matrix& a, vector &b)
{
	int i, j, l;
	matrix s(b.get_size(), b.get_size());
	vector z(b.get_size());
	vector x(b.get_size());
	double temp;

	for (i = 0; i < b.get_size(); i++) {
		for (l = 0, temp = 0; l < i; l++)
			temp += s[l][i]*s[l][i];
		if ((temp = a[i][i]-temp) < 0) {
			cout << "Error: matrix isn't positive.";
			exit(1);
		}
		s[i][i] = sqrt(temp);
		for (j = i+1; j < b.get_size(); j++) {
			for (l = 0, temp = 0; l < i; l++)
				temp += s[l][i]*s[l][j];
			s[i][j] = (a[i][j]-temp)/s[i][i];
		}
	}
	for (i = 0; i < b.get_size(); i++) {
		for (l = 0, temp = 0; l < i; l++)
			temp += s[l][i]*z[l];
		z[i] = (b[i]-temp)/s[i][i];
	}
	for (i = b.get_size()-1; i >= 0; i--) {
		for (l = i+1, temp = 0; l < b.get_size(); l++)
			temp += s[i][l]*x[l];
		x[i] = (z[i]-temp)/s[i][i];
	}
	return x;
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
