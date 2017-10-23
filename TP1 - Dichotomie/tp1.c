#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Prototypes */
double dicho(double x1, double x2, double e, double (*f)(double));
double f1(double x);
double f2(double x);
double complexite(double x1, double x2, double e);

int main()
{
	printf("Racine [1, 2] f1 précision 10^(-5)\n");
	printf("%f\n", dicho(1, 2, pow(10, -5), f1));
	printf("Complexité\n");
	printf("%f\n\n", complexite(1, 2, pow(10, -5)));
	
	printf("Racine [4, 5] f1 précision 10^(-5)\n");
	printf("%f\n", dicho(4, 5, pow(10, -5), f1));
	printf("Complexité\n");
	printf("%f\n\n", complexite(4, 5, pow(10, -5)));
	
	printf("Racine [2, 3] f2 précision 10^(-2)\n");
	printf("%f\n", dicho(2, 3, pow(10, -2), f2));
	printf("Complexité\n");
	printf("%f\n\n", complexite(2, 3, pow(10, -2)));
	
	printf("Racine [2, 3] f2 précision 10^(-3)\n");
	printf("%f\n", dicho(2, 3, pow(10, -3), f2));
	printf("Complexité\n");
	printf("%f\n\n", complexite(2, 3, pow(10, -3)));
	
	printf("Racine [2, 3] f2 précision 10^(-4)\n");
	printf("%f\n", dicho(2, 3, pow(10, -4), f2));
	printf("Complexité\n");
	printf("%f\n\n", complexite(2, 3, pow(10, -4)));
	
	printf("Racine [2, 3] f2 précision 10^(-5)\n");
	printf("%f\n", dicho(2, 3, pow(10, -5), f2));
	printf("Complexité\n");
	printf("%f\n\n", complexite(2, 3, pow(10, -5)));
	
	return 0;
}

double dicho(double x1, double x2, double e, double (*f)(double))
{
	double y1 = (*f)(x1);
	double ym = y1;
	double xm = 0;	
	
	while((x2-x1) > e)
	{
		xm = (x1 + x2) / 2;
		ym = (*f)(xm);
		if ((y1 * ym) < 0)
		{
			x2 = xm; // choix de [x1, xm]
		}
		else
		{
			x1 = xm;
			y1 = ym; // choix de [xm, x2 ]
		}
	}
	return xm;
}

double f1(double x)
{
	return exp(x) - pow(x, 3);
}

double f2(double x)
{
	return pow(x, 2) - 5;
}

double complexite(double x1, double x2, double e)
{
	return (log((x2-x1)/e))/log(2) + 1;
}
