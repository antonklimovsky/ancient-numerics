//--------------------------------------------------
// YASM - Yet Another Simplex Method implementation
// ver. 0.1, May 27, 2001
//
// Author: A. Klimovsky, root@ludus.kharkiv.com
// All rights reserved
//--------------------------------------------------
//--------------------------------------------------
// Simplex method routines
//--------------------------------------------------
#ifndef __GAUSS_H
#define __GAUSS_H

#include <stdio.h>

#include "LinAlg.h"
#include "Log.h"

//////////////////////////////////////////////////////////////////
// Declaration starts here
//////////////////////////////////////////////////////////////////

template <class T>
class SimplexMethod {
private:
    Matrix<T> a;
    Vector<T> b;
    Vector<T> c;
    Vector<T> x;
    Vector<T> delta;
    Vector<T> basisC;
    T f;
    int* basis;
    int* basisVecSet;
    int n;
    int m;
    int k;
    int l;
    bool consistent;
    bool infty;
    bool end;
    Log<T>& log;
    const int bufsize;
    char* buf;

    void init();
    void done();
    void getSolution();
    void calcDelta();
    void checkInfty();
    void chooseBestElement();
    void searchForBasis();
    void includeToBasis();

public:

    SimplexMethod(Matrix<T> a, Vector<T> b, Vector<T> c, Log<T>& log);

    Matrix<T> & getA();

    void setA(Matrix<T> & a);

    Vector<T> & getB();

    void setB(Vector<T> & b);

    Vector<T> & getC();

    void setC(Vector<T> & c);

    Vector <T> & getX();

    bool go();

};

//////////////////////////////////////////////////////////////////
// Implementation starts here
//////////////////////////////////////////////////////////////////

template<class T>
SimplexMethod<T>::SimplexMethod(Matrix<T> a0, Vector<T> b0, Vector<T> c0, Log<T>& log0):
    a(a0),
    b(b0),
    c(c0),
    log(log0),
    bufsize(256)
{
}

template<class T>
Vector<T>& SimplexMethod<T>::getX(){ return x; }

template<class T>
void SimplexMethod<T>::setC(Vector<T>& c){ this->c = c; }

template<class T>
Vector<T>& SimplexMethod<T>::getC(){ return c; }

template<class T>
void SimplexMethod<T>::setB(Vector<T>& b){ this->b = b; }

template<class T>
Vector<T>& SimplexMethod<T>::getB(){ return b; }

template<class T>
void SimplexMethod<T>::setA(Matrix<T>& a){ this->a = a; }

template<class T>
Matrix<T>& SimplexMethod<T>::getA(){ return a; }

template<class T>
void SimplexMethod<T>::searchForBasis()
{
    int i;
    int j;
    int temp;

    for (j = 0; j < m; j++)
        basisVecSet[j] = 0;

    // let's search for the basis columns

    for (i = 0; i < n; i++) {

        basis[i] = -1;

        for (j = 0; j < m && a[j][i] == 0; j++);

        if (j == m || a[j][i] != 1) continue;

        temp = j;

        for (j++; j < m && a[j][i] == 0; j++);

        if (j != m) continue;

        basis[i] = temp;
        basisVecSet[temp] = i+1;
    }
}

template<class T>
void SimplexMethod<T>::init()
{
    int i;
    int j;
    int temp;

    n = a.get_cols();
    m = a.get_rows();
    basis = new int[n];
    basisVecSet = new int[m];
    buf = new char[bufsize];
    consistent = true;

    searchForBasis();

    // now let's check if we've got all of the basis vectors

    for (j = 0; j < m; j++)
        if (!basisVecSet[j]) {
            consistent = false;
            log("[Error] Basis column(s) expected in matrix A\n");
            break;
        }

    // by the way, what about consistency of input data dimensions?

    if (c.get_size() != a.get_cols() || b.get_size() != a.get_rows()) {
        log("[Error] Inconsistent input data dimensions\n");
        consistent = false;
    }

    x = Vector<T>(n);
    delta = Vector<T>(n);
    basisC = Vector<T>(m);
}

template<class T>
void SimplexMethod<T>::done()
{
    delete [] basis;
    delete [] basisVecSet;
    delete [] buf;
}

template<class T>
void SimplexMethod<T>::getSolution()
{
    int i;
    int j;

    log("[Vector c]\n");
    log(c);

    log("[Matrix A]\n");
    log(a);

    log("[Vector b]\n");
    log(b);

    log("[Solution x]\n");

    for (i = 0; i < n; i++)
        x[i] = 0;
    for (j = 0; j < m; j++)
        x[basisVecSet[j]-1] = b[j];

    log(x);
}

template<class T>
void SimplexMethod<T>::calcDelta()
{
    int i;
    int j;
    T temp;

    // let's fill the basisC vector

    for (j = 0; j < m; j++)
        basisC[j] = c[basisVecSet[j]-1];

    log("[Basis c]\n");
    log(basisC);

    // let's calclate the delta's

    for (i = 0; i < n; i++)
        if (basis[i] < 0) {
            temp = c[i];
            for (j = 0; j < m; j++)
                temp -= basisC[j]*a[j][i];
            delta[i] = temp;
        } else
            delta[i] = 0;            

    log("[Delta's]\n");
    log(delta);

    f = 0;
    for (j = 0; j < m; j++)
        f += basisC[j]*b[j];

    log("[Functional f]\n");
    log(f);
    log("\n");
}

template<class T>
void SimplexMethod<T>::checkInfty()
{
    int i;
    int j;

    end = true;
    for (i = 0; i < n; i++)
        if (delta[i] < 0) {
            end = false;
            infty = true;
            for (j = 0; j < m; j++)
                if (a[j][i] > 0) {
                    infty = false;
                    break;
                }
            if (infty) {
                log("Unbounded solution...\n");
                return;
            }
        }
}

template<class T>
void SimplexMethod<T>::chooseBestElement()
{
    int i;
    int j;
    T temp;
    T theta;

    // let's choose the best col

    temp = 0;
    for (i = 0; i < n; i++)
        if (delta[i] < temp) {
            temp = delta[i];
            k = i;
        }

    // let's choose the best row

    // let's skip leading non-positive elements of the best col

    for (j = 0; j < m && a[j][k] <= 0; j++);
    theta = b[j]/a[j][k];
    l = j;
    for (j++; j < m; j++)
        if (a[j][k] > 0 && ((temp = b[j]/a[j][k]) < theta)) {
            theta = temp;
            l = j;
        }
    sprintf(buf, "[Leading element]\n(%i, %i)\n", l+1, k+1);
    log(buf);
}

template<class T>
void SimplexMethod<T>::includeToBasis()
{
    int j;
    int i;
    T temp;

    temp = a[l][k];
    a[l] /= temp;
    b[l] /= temp;
    for (j = 0; j < m; j++)
        if (j != l) {
            temp = a[j][k];
            a[j] -= a[l]*temp;
            b[j] -= b[l]*temp;
        }

    // we've just removed some coord(s) of the solution from basis,
    // but included the other(s), so we have to update our basis[]
    searchForBasis();
}

template<class T>
bool SimplexMethod<T>::go()
{
    int i;

    init();
    i = 0;

    if (consistent)
        do {
            i++;
            sprintf(buf, "[--- Iteration #%i ---]\n", i);
            log(buf);

            getSolution();
            calcDelta();
            checkInfty();
            if (infty)
                break;
            if (end)
                break;
            chooseBestElement();
            includeToBasis();
        } while (true);

    done();

    return !(infty || !consistent);
}

#endif