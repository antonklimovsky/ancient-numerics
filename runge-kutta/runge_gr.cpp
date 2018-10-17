#include <math.h>
#include <stdio.h>
#include <conio.h>
#include <graphics.h>
#include <process.h>
#include <malloc.h>

typedef double (*fp)(double);
int n_points;                  // Number of precalculated points
double x0, x1, y0;
double* points;
double scalex;
double scaley;
double max;
double step;

double f(double x, double y)
{
	return x;
}

double exact_root(double x)
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

void runge_kutt(double** points, double x0, double y0, double x1, double* max,
	 int* n_points, double h)
{
	double x, y;
	int i;

	*n_points = 1+(x1-x0)/h;
	*points = (double *) calloc(*n_points, sizeof(double));
	x = x0;
	**points = y = y0;
	for (i = 1, *max = 0; i < 1+(*n_points); i++) {
		y = make_step(x, y, h);
		x += h;
		(*points)[i] = y;
		if (fabs(y) > *max)
			*max = fabs(y);
	}
}

void draw(fp fun, double scalex, double scaley, double a, double b, double h)
{
	double x;

	moveto(0, (getmaxy()-fun(a)*scaley)/2);
	for (x = x0; x <= x1; x+= h)
		lineto((x-a)*scalex,
				(getmaxy()-fun(x)*scaley)/2);
	getch();
}

void draw0(double scalex, double scaley, double a, double b, double h)
{
	double x;
	int i;

	moveto(0, (getmaxy()-(*points)*scaley)/2);
	for (x = x0, i = 0; x <= x1; x+= h, i++)
		lineto((x-a)*scalex,
				(getmaxy()-points[i]*scaley)/2);
	getch();
}

void main()
{
	int gdriver = DETECT, gmode, errorcode;
	char ch;
	int i;

	printf("\nInput left bound of interval:\nx = ");
	scanf("%lf", &x0);
	printf("y = ");
	scanf("%lf", &y0);
	printf("Input right bound of interval: ");
	scanf("%lf", &x1);

	while (1) {
		printf("\nInput step: ");
		scanf("%lf", &step);

		initgraph(&gdriver, &gmode, "d:\\bc\\bgi");
		runge_kutt(&points, x0, y0, x1, &max, &n_points, step);
		scalex = getmaxx()/(x1-x0);
		scaley = getmaxy()/max;
		draw(exact_root, scalex, scaley, x0, x1, step);
		setcolor(RED);
		draw0(scalex, scaley, x0, x1, step);
		closegraph();
		puts("Continue (y/n)?");
		ch = getch();
		if (ch != 'Y' && ch != 'y')
			break;
		free(points);
	}
}