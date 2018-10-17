#ifndef __MATRIX_H
#define __MATRIX_H

#include <string.h>
#include <iostream.h>

class Matrix {
	double *p;
	long int n;

public:
	Matrix(long int n)
	{   
		this->n = n;
		p = new double[n*n];
	}
	~Matrix()
	{
		delete [] p;
	}
	double& operator() (long int i, long int j)
	{
		return p[i*n+j];
	}
	void operator= (Matrix& lhs)
	{
		memcpy(p, lhs.p, n*n*sizeof(double));
	}
};

#endif