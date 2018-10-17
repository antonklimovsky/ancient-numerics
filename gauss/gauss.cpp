/*
	An implementation of the Gauss algorithm (with selection of the primary element)
	for reversing a given matrix.
*/
#include <iostream.h>
#include <string.h>
#include <math.h>
#include <process.h>
#include "lin_alg.h"
//////////////////////////////////////////////////////////////////
// Algorithm starts here
//////////////////////////////////////////////////////////////////
#define EPS 1.0e-100
int find_largest_in_col(matrix &a, int col)
{
	int size = a.get_rows();
	int i, largest_indx;
	double largest;

	for (i = col+1,
		 largest_indx = col,
		 largest = fabs(a[col][col]); i < size; i++)
		if (fabs(a[i][col]) > largest) {
			largest_indx = i;
			largest = fabs(a[i][col]);
		}
	if (largest <= EPS) {
		cout << "Overflow or singular matrix.\n";
		exit(1);
	}
	return largest_indx;
}
matrix calculate_inverse_matrix(matrix& a)
{
	int size = a.get_rows();
	matrix result(size, size);
	matrix b(size, size << 1);
	vector temp;
	int i, j, best;

// Let's assume B := [A|I], where I is a unit matrix
	for (i = 0; i < size; i++)
		for (j = 0; j < (size << 1); j++)
			if (j < size)
				b[i][j] = a[i][j];
			else
				if ((j-size) == i)
					b[i][j] = 1.0;
				else
					b[i][j] = 0.0;

// Let's triangulate the first block of B
	for (i = 0; i < size-1; i++)
	{
		best = find_largest_in_col(b, i);
		b.swap_rows(i, best);
		temp = (1.0/b[i][i])*b[i];
		for (j = i+1; j < size; j++)
			if (fabs(b[j][i]) > EPS)
				b[j] -= b[j][i]*temp;

	}

// Let's check if the last row is a zero vector
	if (fabs(b[i][i]) <= EPS) {
		cout << "Overflow or singular matrix.\n";
		exit(1);
	}

// Let's transform B to the form [I|A^(-1)]
	for (i = size-1; i >= 0; i--) {
		b[i] *= 1.0/b[i][i];
		for (j = i-1; j >= 0; j--)
			b[j] -= b[j][i]*b[i];
	}

// Now B = [I|A^(-1)], so RESULT := A^(-1)
	for (i = 0; i < size; i++)
		for (j = 0; j < size; j++)
			result[i][j] = b[i][size+j];
	return result;
}
matrix get_matrix()
{
	int size;

	cout << "Input matrix dimension: ";
	cin >> size;
	matrix a(size, size);
	cout << "Input matrix:\n";
	cin >> a;
	return a;
}
/*
//////////////////////////////////////////////////////////////////
// Following function returns A(i,j):=1/(i+j+1), where i,j=0..size
//////////////////////////////////////////////////////////////////
matrix get_matrix()
{
	int size;
	int i, j;

	cout << "Input matrix dimension: ";
	cin >> size;
	matrix a(size, size);
	for (i = 0; i < size; i++)
		for (j = 0; j < size; j++)
			a[i][j] = 1.0/(i+j+1.0);
	cout << a;
	return a;
}*/
//////////////////////////////////////////////////////////////////
// Main body of the program
//////////////////////////////////////////////////////////////////
void main()
{
	matrix a, b;

	a = get_matrix();
	b = calculate_inverse_matrix(a);
	cout << "Inverse matrix:\n" << b;
	cout << "Product:\n" << a*b;
}
