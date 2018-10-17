////////////////////////////////////////////////////////////////////////////////////
// A a simple abstraction to work with EXACT fractions of numbers of a given type T
////////////////////////////////////////////////////////////////////////////////////

#ifndef __FRACT_H
#define __FRACT_H

#include <math.h>

template<class T> class Fraction {
private:

  T nom;
  T denom;

  T nod(T n, T d);

public:
    Fraction():
        nom(0),
        denom(1)
    {}

    Fraction(T n, T d = 1) {
        SetFraction(n,d);
    }

    void setFraction(T n, T d = 1) {
        if (n == 0) { nom = 0; return; }
        T z=nod(n,d);
        nom=(d/abs(d))*(n/z);
        denom=abs(d)/z;
    }

    T getNom() {
        return nom;
    }

    T getDeNom() {
        return denom;
    }

    double toDouble() {
        return (*this);
    }

    Fraction<T> operator* (Fraction<T> a) {
        return Fraction<T>(nom*a.GetNom(),denom*a.GetDeNom());
    }

    Fraction<T> operator/ (Fraction<T> a) {
        return Fraction<T>(nom*a.GetDeNom(),denom*a.GetNom());
    }

    Fraction<T> operator+ (Fraction<T> a) {
        return Fraction<T>(nom*a.GetDeNom()+denom*a.GetNom(),denom*a.GetDeNom());
    }

    Fraction<T> operator- (Fraction<T> a) {
        return Fraction<T>(nom*a.GetDeNom()-denom*a.GetNom(),denom*a.GetDeNom());
    }

    operator double () {
        return ((double) nom)/denom;
    }

}

template<class T> T Fraction<T>::nod(T n, T d) {
    if (d == 0) return 1; // ???
    if (n == 0) return 0;
    T c; n=abs(n); d=abs(d);
    if (d>n) {c=n; n=d; d=c;};
    while ((c=n%d) != 0) {n=d; d=c;}
    return d;
}

typedef Fraction<int> FrInt;

#endif
