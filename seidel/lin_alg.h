#ifndef __lin_alg__
#define __lin_alg__
#include <iostream.h>
///////////////////////////////////////////////////////////////////
//	Vector class definition
///////////////////////////////////////////////////////////////////
class vector {
	int n;
	double *data;

public:
	vector();
	vector(int size);
	vector(vector&);
	~vector();
	int get_size();
	friend ostream& operator<<(ostream&, vector&);
	friend istream& operator>>(istream&, vector&);
	friend vector operator*(double, vector&);
	double operator*(vector&);
	vector operator+(vector&);
	vector operator-(vector&);
	vector operator=(vector&);
	vector operator+=(vector&);
	vector operator-=(vector&);
	vector operator*=(double);
	vector operator-();
	vector& operator+();
	double& operator[](int i);
};
///////////////////////////////////////////////////////////////////
//	Matrix class definition
///////////////////////////////////////////////////////////////////
class matrix {
	int r, c;
	vector** data;

public:
	matrix();
	matrix(matrix&);
	matrix(int rows, int cols);
	~matrix();
	void swap_rows(int, int);
	int get_rows();
	int get_cols();
	friend ostream& operator<<(ostream&, matrix&);
	friend istream& operator>>(istream&, matrix&);
	vector& operator[](int);
	matrix& operator=(matrix&);
	matrix operator*(matrix&);
	vector operator*(vector&);
};
double abs(vector);
#endif