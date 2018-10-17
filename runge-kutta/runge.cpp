/*
*/
#include <iostream.h>
#include <math.h>

double f(double x, double y)
{
	return x;
}

double f0(double x)
{
	return x*x/2;
}

double make_step(double x, double y, double h)
{
	double k1, k2, k3, k4;

	k1 = f(x, y);
	k2 = f(x+h/2, y+h*k1/2);
	k3 = f(x+h/2, y+h*k2/2);
	k4 = f(x+h, y+h*k3);
	return y+h*(k1+2*k2+2*k3+k4)/6;
}

void main()
{
	double x;
	double y;
	double h;
	int n, i;

	cout << "Input starting point:\nx = ";
	cin >> x;
	cout << "y = ";
	cin >> y;
	cout << "Step length: ";
	cin >> h;
	cout << "Number of steps: ";
	cin >> n;

	for (i = 0; i < n; i++) {
		y = make_step(x, y, h);
		x += h;
		cout << "f(" << x << ") = " << y << ", " << f0(x) << "\n";
	}
}