#ifndef __PROGONKA_H
#define __PROGONKA_H

class LinearSystem { // with 3-diagonal matrix
public:
	long int n;
// Coefficients
	double *a;
	double *b;
	double *c;
	double *f;
	double k1, k2;
	double n1, n2;
// Methods
	LinearSystem(int n);
	~LinearSystem();
	double* solve();
};
#endif
