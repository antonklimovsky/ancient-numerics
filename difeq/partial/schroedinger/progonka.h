#ifndef __PROGONKA_H
#define __PROGONKA_H

namespace Progonka {
	template<class T> class LinearSystem {
	public:
		long int n;

// Coefficients of the equation

		T *a;
		T *b;
		T *c;
		T *f;
		T* x;
		T k1, k2;
		T n1, n2;

		LinearSystem<T>(long int n0)
		{
			n = n0;
			a = new T[n-1];
			b = new T[n-1];
			c = new T[n-1];
			f = new T[n-1];
			x = new T[n+1];
		}
		~LinearSystem()
		{
			delete [] a;
			delete [] b;
			delete [] c;
			delete [] f;
			delete [] x;
		}
		T* solve()
		{
			T* alpha;
			T* beta;
			long int i;

			if(!(alpha = new T[n])) { cout << "Not enough memory"; exit(1); };
			if(!(beta = new T[n])) { cout << "Not enough memory"; exit(1); };

			*alpha = k1;
			*beta = n1;
			for (i = 0; i < n-1; i++) {
				alpha[i+1] = -c[i]/(b[i]+a[i]*alpha[i]);
				beta[i+1] = (f[i]-a[i]*beta[i])/(b[i]+a[i]*alpha[i]);
			}

			x[n] = (k2*beta[n-1]+n2)/(complex<double>(1, 0)-k2*alpha[n-1]);
			for (i = n-1; i >= 0; i--)
				x[i] = alpha[i]*x[i+1]+beta[i];
			delete [] alpha;
			delete [] beta;
			return x;
		}
	};
};

#endif
