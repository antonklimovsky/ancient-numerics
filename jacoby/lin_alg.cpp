#include <math.h>
#include "lin_alg.h"
///////////////////////////////////////////////////////////////////
//	Vector class implementation
///////////////////////////////////////////////////////////////////
vector::vector()
{
	n = 0;
}
inline vector::vector(int size)
{
	data = new double[size];
	n = size;
}
vector::vector(vector& vec)
{
	data = new double[n = vec.n];
	memcpy(data, vec.data, n*sizeof(double));
}
inline vector::~vector()
{
	if (n != 0)
		delete data;
}
inline double& vector::operator[](int i)
{
	return data[i];
}
inline vector& vector::operator+()
{
	return *this;
}
int vector::get_size()
{
	return n;
}
ostream& operator<<(ostream& output, vector& vec)
{
	int i;

	for (i = 0; i < vec.n; i++)
		output << vec.data[i] << " ";
	return output;
}
istream& operator>>(istream& input, vector& vec)
{
	int i;

	for (i = 0; i < vec.n; i++)
		input >> vec.data[i];
	return input;
}
vector operator*(double scalar, vector& vec)
{
	int i;
	vector temp(vec);

	for (i = 0; i < vec.n; i++)
		temp[i] *= scalar;
	return temp;
}
double vector::operator*(vector& vec)
{
	int i;
	double product = 0.0;

	for (i = 0; i < vec.n; i++)
		product += data[i]*vec[i];
	return product;
}
vector vector::operator-()
{
	int i;
	vector temp(*this);

	for (i = 0; i < n; i++)
		temp[i] = -data[i];
	return temp;
}
vector vector::operator+(vector& vec)
{
	int i;
	vector temp(*this);

	for (i = 0; i < n; i++)
		temp[i] += vec.data[i];

	return temp;
}
vector vector::operator-(vector& vec)
{
	int i;
	vector temp(*this);

	for (i = 0; i < n; i++)
		temp[i] -= vec.data[i];

	return temp;
}
vector vector::operator=(vector& vec)
{
	int i;

	if (n != 0)
		delete data;
	data = new double[n = vec.n];
	data = (double *) memcpy(data, vec.data, vec.n*sizeof(double));
	return *this;
}
vector vector::operator+=(vector& vec)
{
	int i;

	if (n != vec.n) return *this;
	for (i = 0; i < n; i++)
		data[i] += vec[i];
	return *this;
}
vector vector::operator-=(vector& vec)
{
	int i;

	if (n != vec.n) return *this;
	for (i = 0; i < n; i++)
		data[i] -= vec[i];
	return *this;
}
vector vector::operator*=(double scalar)
{
	int i;

	for (i = 0; i < n; i++)
		data[i] *= scalar;
	return *this;
}
///////////////////////////////////////////////////////////////////
//	Function abs
///////////////////////////////////////////////////////////////////
double abs(vector a)
{
	double max;
	int i;

	for (i = 0, max = 0; i < a.get_size(); i++)
		if (max < fabs(a[i]))
			max = fabs(a[i]);
	return max;
}
///////////////////////////////////////////////////////////////////
//	Matrix class implementation
///////////////////////////////////////////////////////////////////
matrix::matrix()
{
	r = 0;
	c = 0;
}
matrix::matrix(matrix& mat)
{
	int i;

	r = mat.r;
	c = mat.c;
	data = new vector*[r];
	for (i = 0; i < r; i++) {
		data[i] = new vector(c);
		*data[i] = mat[i];
	}
}
matrix::matrix(int rows, int cols)
{
	int i;

	r = rows;
	c = cols;
	data = new vector*[rows];
	for (i = 0; i < rows; i++)
		data[i] = new vector(cols);
}
matrix::~matrix()
{
	int i;

	for (i = 0; i < r; i++)
		delete data[i];
	delete [] data;
}
inline void matrix::swap_rows(int r1, int r2)
{
	vector *temp;

	temp = data[r1];
	data[r1] = data[r2];
	data[r2] = temp;
}
int matrix::get_rows()
{
	return r;
}
int matrix::get_cols()
{
	return c;
}
ostream& operator<<(ostream& output, matrix& mat)
{
	int i;

	for (i = 0; i < mat.r; i++)
		output << *(mat.data[i]) << "\n";
	return output;
}
istream& operator>>(istream& input, matrix& mat)
{
	int i;

	for (i = 0; i < mat.r; i++)
		input >> *(mat.data[i]);
	return input;
}
vector& matrix::operator[](int row)
{
	return *(data[row]);
}
matrix& matrix::operator=(matrix& mat)
{
	int i;

	if (r != 0) {
		for (i = 0; i < r; i++)
			delete data[i];
		delete [] data;
	}
	r = mat.r;
	c = mat.c;
	data = new vector*[r];
	for (i = 0; i < r; i++) {
		data[i] = new vector(c);
		*data[i] = mat[i];
	}
	return *this;
}
matrix matrix::operator*(matrix& mat)
{
	int i, j, k;
	matrix product(r, mat.c);
	double temp, t1;

	if (c == mat.r)
		for (i = 0; i < r; i++)
			for (j = 0; j < mat.c; j++) {
				temp = 0;
				for (k = 0; k < c; k++)
					temp += (*this)[i][k]*mat[k][j];
				product[i][j] = temp;
			}
	return product;
}
vector matrix::operator*(vector& vec)
{
	int i, j;
	vector product(vec.get_size());

	for (i = 0; i < vec.get_size(); i++) {
		product[i] = 0;
		for (j = 0; j < vec.get_size(); j++)
			product[i] += (*this)[i][j]*vec[j];
	}
	return product;
}
matrix matrix::transpose()
{
	int i, j;
	matrix temp(r, r);

	for (i = 0; i < r; i++)
		for (j = 0; j < r; j++)
			temp[i][j] = (*this)[j][i];
	return temp;
}