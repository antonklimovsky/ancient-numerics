//--------------------------------------------------
// YASM - Yet Another Simplex Method implementation
// ver. 0.1, May 27, 2001
//
// Author: A.Klimovsky, root@ludus.kharkiv.com
// Special thanks to M.Zinchenko, mxz@inbox.ru
// All rights reserved
//--------------------------------------------------
//-----------------------------------------
// The rational number abstraction routines
//-----------------------------------------
#ifndef __RATIONAL_H
#define __RATIONAL_H

#include <math.h>
#include <stdio.h>

#include "miracl/big.h"
#include "miracl/flash.h"

template<class T> class Rational {
private:
  T nom;
  T denom;

  T nod(T n, T d);

public:
    Rational():
        nom(0),
        denom(1)
    {}

    Rational(T n, T d = 1)
    {
        setRational(n,d);
    }

    void setRational(T n, T d = 1) {
        T z;
        if (n == 0) {
            nom = 0;
            denom = 1;
            return;
        }
        z = nod(n,d);
        nom = (d/abs(d))*(n/z);
        denom=abs(d)/z;
    }

    T getNom() {
        return nom;
    }

    T getDenom() {
        return denom;
    }

    double toDouble() {
        return (*this);
    }

    Rational<T> operator* (Rational<T> a) {
        Rational<T> temp(nom, a.denom);
        Rational<T> temp1(a.nom, denom);

        return Rational<T>(temp.nom*temp1.nom,temp.denom*temp1.denom);
    }

    Rational<T> operator/ (Rational<T> a) {
        Rational<T> temp(nom, a.nom);
        Rational<T> temp1(a.denom, denom);

        return Rational<T>(temp.nom*temp1.nom,temp.denom*temp1.denom);
    }

    Rational<T> operator+ (Rational<T> a) {
        T z = nod(denom, a.denom);
        T aa = a.denom/z;
        T bb = denom/z;
        return Rational<T>(nom*aa+a.nom*bb,aa*denom);
    }

    Rational<T> operator- (Rational<T> a) {
        T z = nod(denom, a.denom);
        T aa = a.denom/z;
        T bb = denom/z;
        return Rational<T>(nom*aa-a.nom*bb,aa*denom);
    }

    Rational<T> operator= (Rational<T> a) {
        nom = a.getNom();
        denom = a.getDenom();
        return *this;
    }

    Rational<T> operator= (T a) {
        nom = a;
        denom = 1;
        return *this;
    }

    void operator/= (Rational<T> a) {
        Rational<T> temp(nom, a.nom);
        Rational<T> temp1(a.denom, denom);

        setRational(temp.nom*temp1.nom,temp.denom*temp1.denom);
    }

    void operator*= (Rational<T> a) {
        Rational<T> temp(nom, a.denom);
        Rational<T> temp1(a.nom, denom);

        setRational(temp.nom*temp1.nom,temp.denom*temp1.denom);
    }

    void operator+= (Rational<T> a) {
        T z = nod(denom, a.denom);
        T aa = a.denom/z;
        T bb = denom/z;
        setRational(nom*aa+a.nom*bb,aa*denom);
    }

    void operator-= (Rational<T> a) {
        T z = nod(denom, a.denom);
        T aa = a.denom/z;
		T bb = denom/z;
		setRational(nom*aa-a.nom*bb,aa*denom);
	}

	/*operator double () {
		return todouble(Flash(nom)/Flash(denom));
	}*/

	bool operator< (Rational<T> a) {
		Rational<T> temp =  *this - a;
		return temp.nom < 0;
	}

	bool operator< (long a) {
		return (nom-denom*a) < 0;
	}

	bool operator<= (Rational<T> a) {
		Rational<T> temp =  *this - a;
		return temp.nom <= 0;
	}

	bool operator<= (long a) {
		return (nom-denom*a) <= 0;
	}

	bool operator== (Rational<T> a) {
		Rational<T> temp =  *this - a;
		return temp.nom == 0;
	}

	bool operator== (long a) {
		return (nom-denom*a) == 0;
	}

	bool operator!= (Rational<T> a) {
		Rational<T> temp =  *this - a;
		return temp.nom != 0;
	}

	bool operator!= (long a) {
		return (nom-a) != 0;
	}

	bool operator> (Rational<T> a) {
		Rational<T> temp =  *this - a;
		return temp.nom > 0;
	}

	bool operator> (long a) {
		return (nom-denom*a) > 0;
	}

	bool operator>= (Rational<T> a) {
		Rational<T> temp =  *this - a;
		return temp.nom >= 0;
	}

	bool operator>= (long a) {
		return (nom-denom*a) >= 0;
	}
};

template<class T> T Rational<T>::nod(T n, T d)
{
    if (d == 0)
        return 1; // ???
    if (n == 0)
        return 1;
    T c;
    n=abs(n);
    d=abs(d);
    if (d > n) {
        c=n;
        n=d;
        d=c;
    }
    while ((c=n%d) != 0) {
        n=d;
        d=c;
    }
    return d;
}

template<class T> ostream& operator<<(ostream& output, Rational<T> r)
{
    const int BUFSIZE = 256;
    char buf[BUFSIZE];

    output.width(7);
    if (r.getNom() == 0) {
        output << "0";
    }
    else if (r.getDenom() == 1) {
        output << r.getNom();
    }
    else {
        sprintf(buf, "%i/%i", r.getNom(), r.getDenom());
        output << buf;
    }
    return output;
}

typedef Rational<Big> FractInt;

#endif