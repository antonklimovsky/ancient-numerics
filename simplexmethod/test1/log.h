//--------------------------------------------------
// YASM - Yet Another Simplex Method implementation
// ver. 0.1, May 27, 2001
//
// Author: A.Klimovsky, root@ludus.kharkiv.com
// All rights reserved
//--------------------------------------------------
//--------------------------------------------------
// The extremely simple logger
//--------------------------------------------------
#ifndef __LOG_H
#define __LOG_H

#include "Rational.h"
#include "LinAlg.h"

template<class T> class Log {
public:
    Log() {}
    void operator() (char* s)
    {
        cout << s;
    }
    void operator() (T r)
    {
        cout << r << " ";
    }
    void operator() (Vector<T> v)
    {
        cout << v << "\n";
    }
    void operator() (Matrix<T> m)
    {
        cout << m;
    }
};

#endif