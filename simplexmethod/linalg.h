//--------------------------------------------------
// YASM - Yet Another Simplex Method implementation
// ver. 0.1, May 27, 2001
//
// Author: A. Klimovsky, root@ludus.kharkiv.com
// All rights reserved
//--------------------------------------------------
//--------------------------------------------------
// Some linear algebra routines
//--------------------------------------------------
#ifndef __lin_alg__
#define __lin_alg__

#include <iostream.h>
#include <math.h>

///////////////////////////////////////////////////////////////////
//  Vector class definition
///////////////////////////////////////////////////////////////////
template<class T> class Vector {
    int n;
    T *data;

public:
    Vector();
    Vector(int size);
    Vector(int size, T*);
    Vector(Vector<T>&);
    ~Vector();
    int get_size();
    friend Vector<T> operator* <T>(T, Vector<T>);
    Vector operator*(T);
    Vector operator/(T);
    T operator*(Vector<T>);
    Vector operator+(Vector<T>);
    Vector operator-(Vector<T>);
    Vector operator=(Vector<T>);
    Vector operator+=(Vector<T>);
    Vector operator-=(Vector<T>);
    Vector operator*=(T);
    Vector operator/=(T);
    Vector operator-();
    Vector& operator+();
    T& operator[](int i);
};
///////////////////////////////////////////////////////////////////
//  Matrix class definition
///////////////////////////////////////////////////////////////////
template<class T> class Matrix {
    int r, c;
    Vector<T>** data;

public:
    Matrix();
    Matrix(Matrix<T>&);
    Matrix(int rows, int cols);
    Matrix(int rows, int cols, T* data);
    ~Matrix();
    void swap_rows(int, int);
    int get_rows();
    int get_cols();
    Vector<T>& operator[](int);
    Matrix<T>& operator=(Matrix<T>);
    Matrix<T> operator-();
    Matrix<T> operator+(Matrix<T>);
    Matrix<T> operator-(Matrix<T>);
    Matrix<T> operator*(Matrix<T>);
    Matrix<T> operator/(T);
    friend Matrix<T> operator* <T>(T, Matrix<T>);
    Vector<T> operator*(Vector<T>);
    Matrix<T> transpose();
};

template<class T> T abs(Vector<T>);

///////////////////////////////////////////////////////////////////
//  Vector class implementation
///////////////////////////////////////////////////////////////////
template<class T> Vector<T>::Vector()
{
    n = 0;
}
template<class T> Vector<T>::Vector(int size)
{
    data = new T[size];
    n = size;
}
template<class T> Vector<T>::Vector(int size, T* data)
{
    this->data = new T[size];
    n = size;
    for (int i = 0; i < n; i++)
        this->data[i] = data[i];
}
template<class T> Vector<T>::Vector(Vector<T>& vec)
{
    data = new T[n = vec.n];
    for (int i = 0; i < n; i++)
        data[i] = vec.data[i];
}
template<class T> Vector<T>::~Vector()
{
    if (n != 0)
        delete [] data;
}
template<class T> inline T& Vector<T>::operator[](int i)
{
    return data[i];
}
template<class T> inline Vector<T>& Vector<T>::operator+()
{
    return *this;
}
template<class T> inline int Vector<T>::get_size()
{
    return n;
}
template<class T> ostream& operator<<(ostream& output, Vector<T>& vec)
{
    int i;

    for (i = 0; i < vec.get_size(); i++)
        output << vec[i] << " ";
    return output;
}
template<class T> istream& operator>>(istream& input, Vector<T>& vec)
{
    int i;

    for (i = 0; i < vec.get_size(); i++)
        input >> vec[i];
    return input;
}
template<class T> Vector<T> operator*(T scalar, Vector<T> vec)
{
    int i;
    Vector<T> temp(vec);

    for (i = 0; i < vec.n; i++)
        temp[i] *= scalar;
    return temp;
}
template<class T> Vector<T> Vector<T>::operator*(T scalar)
{
    int i;
    Vector<T> temp(*this);

    for (i = 0; i < n; i++)
        temp[i] *= scalar;
    return temp;
}
template<class T> Vector<T> Vector<T>::operator/(T scalar)
{
    int i;
    Vector<T> temp(*this);

    for (i = 0; i < n; i++)
        temp[i] /= scalar;
    return temp;
}
template<class T> T Vector<T>::operator*(Vector<T> vec)
{
    int i;
    T product = 0.0;

    for (i = 0; i < vec.n; i++)
        product += data[i]*vec[i];
    return product;
}
template<class T> Vector<T> Vector<T>::operator-()
{
    int i;
    Vector<T> temp(*this);

    for (i = 0; i < n; i++)
        temp[i] = -data[i];
    return temp;
}
template<class T> Vector<T> Vector<T>::operator+(Vector<T> vec)
{
    int i;
    Vector<T> temp(*this);

    for (i = 0; i < n; i++)
        temp[i] += vec.data[i];

    return temp;
}
template<class T> Vector<T> Vector<T>::operator-(Vector<T> vec)
{
    int i;
    Vector<T> temp(*this);

    for (i = 0; i < n; i++)
        temp[i] -= vec.data[i];

    return temp;
}
template<class T> Vector<T> Vector<T>::operator=(Vector<T> vec)
{
    if (n != 0)
        delete [] data;
    data = new T[n = vec.n];
    for (int i = 0; i < n; i++)
        data[i] = vec.data[i];   
    return *this;
}
template<class T> Vector<T> Vector<T>::operator+=(Vector<T> vec)
{
    int i;

    if (n != vec.n) return *this;
    for (i = 0; i < n; i++)
        data[i] += vec[i];
    return *this;
}
template<class T> Vector<T> Vector<T>::operator-=(Vector<T> vec)
{
    int i;

    if (n != vec.n) return *this;
    for (i = 0; i < n; i++)
        data[i] -= vec[i];
    return *this;
}
template<class T> Vector<T> Vector<T>::operator*=(T scalar)
{
    int i;

    for (i = 0; i < n; i++)
        data[i] *= scalar;
    return *this;
}
template<class T> Vector<T> Vector<T>::operator/=(T scalar)
{
    int i;

    for (i = 0; i < n; i++)
        data[i] /= scalar;
    return *this;
}
///////////////////////////////////////////////////////////////////
//  Function abs
///////////////////////////////////////////////////////////////////
template<class T> T abs(Vector<T> a)
{
    T max;
    int i;

    for (i = 0, max = 0; i < a.get_size(); i++)
        if (max < fabs(a[i]))
            max = fabs(a[i]);
    return max;
}
///////////////////////////////////////////////////////////////////
//  Matrix class implementation
///////////////////////////////////////////////////////////////////
template<class T> Matrix<T>::Matrix()
{
    r = 0;
    c = 0;
}
template<class T> Matrix<T>::Matrix(Matrix<T>& mat)
{
    int i;

    r = mat.r;
    c = mat.c;
    data = new Vector<T>*[r];
    for (i = 0; i < r; i++) {
        data[i] = new Vector<T>(c);
        *data[i] = mat[i];
    }
}
template<class T> Matrix<T>::Matrix(int rows, int cols)
{
    int i;

    r = rows;
    c = cols;
    data = new Vector<T>*[rows];
    for (i = 0; i < rows; i++)
        data[i] = new Vector<T>(cols);
}
template<class T> Matrix<T>::Matrix(int rows, int cols, T* d)
{
    int i, j;
    Vector<T> temp;

    r = rows;
    c = cols;
    data = new Vector<T>*[rows];
    for (i = 0; i < rows; i++) {
        data[i] = new Vector<T>(cols);
        for (j = 0; j < cols; j++)
            (*data[i])[j] = *(d+i*cols+j);
    }
}
template<class T> Matrix<T>::~Matrix()
{
    int i;

    for (i = 0; i < r; i++)
        delete data[i];
    delete [] data;
}
template<class T> void Matrix<T>::swap_rows(int r1, int r2)
{
    Vector<T> *temp;

    temp = data[r1];
    data[r1] = data[r2];
    data[r2] = temp;
}
template<class T> inline int Matrix<T>::get_rows()
{
    return r;
}
template<class T> inline int Matrix<T>::get_cols()
{
    return c;
}
template<class T> ostream& operator<<(ostream& output, Matrix<T>& mat)
{
    int i;

    for (i = 0; i < mat.get_rows(); i++)
        output << mat[i] << "\n";
    return output;
}
template<class T> istream& operator>>(istream& input, Matrix<T>& mat)
{
    int i;

    for (i = 0; i < mat.get_rows(); i++)
        input >> mat[i];
    return input;
}
template<class T> inline Vector<T>& Matrix<T>::operator[](int row)
{
    return *(data[row]);
}
template<class T> Matrix<T>& Matrix<T>::operator=(Matrix mat)
{
    int i;

    if (r != 0) {
        for (i = 0; i < r; i++)
            delete data[i];
        delete [] data;
    }
    r = mat.r;
    c = mat.c;
    data = new Vector<T>*[r];
    for (i = 0; i < r; i++) {
        data[i] = new Vector<T>(c);
        *data[i] = mat[i];
    }
    return *this;
}
template<class T> Matrix<T> Matrix<T>::operator-()
{
    int i, j;
    Matrix<T> temp(r, c);

    for (i = 0; i < r; i++)
        for (j = 0; j < c; j++)
            temp[i][j] = -(*this)[i][j];
    return temp;
}
template<class T> Matrix<T> Matrix<T>::operator+(Matrix rhs)
{
    int i, j;
    Matrix<T> temp(r, c);

    for (i = 0; i < r; i++)
        temp[i] = (*this)[i]+rhs[i];
    return temp;
}
template<class T> Matrix<T> Matrix<T>::operator-(Matrix<T> rhs)
{
    int i, j;
    Matrix<T> temp(r, c);

    for (i = 0; i < r; i++)
        temp[i] = (*this)[i]-rhs[i];
    return temp;
}
template<class T> Matrix<T> operator*(T d, Matrix<T> mat)
{
    int i, j;
    Matrix<T> temp(mat.r, mat.c);

    for (i = 0; i < mat.r; i++)
        temp[i] = d*mat[i];
    return temp;
}
template<class T> Matrix<T> Matrix<T>::operator/(T d)
{
    int i, j;
    Matrix<T> temp(*this);

    for (i = 0; i < this->r; i++)
        temp[i] = temp[i]/d;
    return temp;
}
template<class T> Matrix<T> Matrix<T>::operator*(Matrix<T> mat)
{
    int i, j, k;
    Matrix<T> product(r, mat.c);
    T temp, t1;

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
template<class T> Vector<T> Matrix<T>::operator*(Vector<T> vec)
{
    int i, j;
    Vector<T> product(vec.get_size());

    for (i = 0; i < vec.get_size(); i++) {
        product[i] = 0;
        for (j = 0; j < vec.get_size(); j++)
            product[i] += (*this)[i][j]*vec[j];
    }
    return product;
}
template<class T> Matrix<T> Matrix<T>::transpose()
{
    int i, j;
    Matrix<T> temp(c, r);

    for (i = 0; i < c; i++)
        for (j = 0; j < r; j++)
            temp[i][j] = (*this)[j][i];
    return temp;
}

#endif