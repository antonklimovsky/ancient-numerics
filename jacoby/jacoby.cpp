/*
	Implementation of the Jacoby rotation algorithm for finding the spectrum of a
	symmetric matrix.
*/
#include <iostream.h>
#include <string.h>
#include <math.h>
#include <process.h>
#include "lin_alg.h"
//////////////////////////////////////////////////////////////////
// Algorithm starts here
//////////////////////////////////////////////////////////////////
#define MAX_ITERATIONS 100
#define EPS 1.0e-5
double find_largest(matrix& a, int& best_i, int& best_j)
{
	int i, j;
	double max = 0;

	best_i = 0;
	best_j = 0;
	for (i = 0; i < a.get_rows(); i++)
		for (j = i+1; j < a.get_rows(); j++)
			if (fabs(a[i][j]) > max) {
				best_i = i;
				best_j = j;
				max = fabs(a[i][j]);
			}
	return max;
}
double inline sqr(double x)
{
	return x*x;
}
double inline sign(double x)
{
	return (x >= 0)? 1 : -1;
}
void calculate_csd(matrix& a, int i, int j, double& c, double& s, double& d)
{
	d = sqrt(sqr(a[i][i]-a[j][j])+4*sqr(a[i][j]));
	c = sqrt((1+fabs(a[i][i]-a[j][j])/d)/2);
	s = sign(a[i][j]*(a[i][i]-a[j][j]))*sqrt((1-fabs(a[i][i]-a[j][j])/d)/2);
}
void make_id(matrix& a)
{
	int i, j;

	for (i = 0; i < a.get_rows(); i++)
		for (j = 0; j < a.get_rows(); j++)
			if (i == j)
				a[i][j] = 1;
			else
				a[i][j] = 0;
}
void produce_u(matrix& u, int i, int j, double c, double s)
{
	make_id(u);
	u[i][i] = c;
	u[j][j] = c;
	u[i][j] = -s;
	u[j][i] = s;
}
void iterate(matrix& a, matrix& temp, int i, int j, double c, double s,
	double d)
{
	int k, l;

	for (k = 0; k < a.get_rows(); k++)
		if (k  != i && k != j)
			for (l = 0; l < a.get_rows(); l++)
				if (l  != i && l != j)
					temp[k][l] = a[k][l];
	for (k = 0; k < a.get_rows(); k++)
		if (k  != i && k != j) {
			temp[k][i] = temp[i][k] = c*a[k][i]+s*a[k][j];
			temp[k][j] = temp[j][k] = -s*a[k][i]+c*a[k][j];
		}
   /*	temp[i][i] = (a[i][i]-a[j][j])/2+sign(a[i][i]-a[j][j])*d/2;
	temp[j][j] = (a[i][i]+a[j][j])/2-sign(a[i][i]-a[j][j])*d/2;*/
	temp[i][i] = sqr(c)*a[i][i]+2*c*s*a[i][j]+sqr(s)*a[j][j];
	temp[j][j] = sqr(s)*a[i][i]-2*c*s*a[i][j]+sqr(c)*a[j][j];
	temp[i][j] = temp[j][i] = 0;
}
vector find_spectrum(matrix& a, matrix &basis_)
{
	int i, j;
	int best_i, best_j;
	double max;
	vector spectrum(a.get_rows());
	matrix basis(a.get_rows(), a.get_rows());
	matrix u(a.get_rows(), a.get_rows());
	matrix c(a.get_rows(), a.get_rows());
	matrix temp(a.get_rows(), a.get_rows());
	double cosine, sine, d;


	make_id(basis);
	c = a;
	for (i= 0; i < MAX_ITERATIONS; i++) {
		if (find_largest(c, best_i, best_j) < EPS) {
			cout << i << " iterations have been done.\n";
			basis_ = basis;
			for (j = 0; j < a.get_rows(); j++)
				spectrum[j] = c[j][j];
			return spectrum;
		}
		calculate_csd(c, best_i, best_j, cosine, sine, d);
		produce_u(u, best_i, best_j, cosine, sine);
		basis = basis*u;
		iterate(c, temp, best_i, best_j, cosine, sine, d);
		c = temp;
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
	vector spec;
	matrix basis;
	int n;

	cout << "Input dimension of the matrix: ";
	cin >> n;

	cout << "Input matrix:\n";
	a = matrix(n, n);
	cin >> a;

	spec = find_spectrum(a, basis);
	cout << "Spectrum:\n" << spec << "\n";
	cout << "Basis:\n" << basis;
	cout << "Test:\n" << basis.transpose()*a*basis << "\n";
}
