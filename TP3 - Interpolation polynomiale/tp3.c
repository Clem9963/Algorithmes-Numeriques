#include <stdio.h>
#include <stdlib.h>

// PROTOTYPES

double lagrange(double *x, double *y, int n, double a);
double newton(double *x, double *y, int n, double a);

int main()
{
	double x[20] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38};
	double y[20] = {0.99987, 0.99997, 1, 0.99997, 0.99988, 0.99973, 0.99953, 0.99927, 0.99897, 0.99846,
					0.99805, 0.99751, 0.99705, 0.99650, 0.99664, 0.99533, 0.99472, 0.99472, 0.99333, 0.99326};

	printf("%f\n", lagrange(x, y, 20, 0));
	return 0;
}

double lagrange(double *x, double *y, int n, double a)
{
	int i = 0;
	int j = 0;
	double li = 1;
	double somme = 0;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (j != i)
			{
				li = li * ((a-x[j])/(x[i]-x[j]));
			}
		}
		somme = somme + y[i] * li;
		li = 1;
	}
	return somme;
}

double newton(double *x, double *y, int n, double a)
{
	int a = 0;
	int b = 0;
	int ordre = 0;
	double* differences_divisees[n-1] = {NULL};
	for (a = 0; a < n-1; a++)
	{
		differences_divisees[a] = calloc(n-1-a, sizeof(double));
	}
	for (a = 0; a < n-1; a++)
	{
		for (b = 0; b < n-1-a; b++)
		{
			ordre = a + 1;
			differences_divisees[a][b] =  /
		}
	}
}

double trouve_a(double **differences_divisees, double *a, int n)
{
	int i = 0;
	int j = 0;
	for (i = 0; i < n; i++)
	{
		a[i] = differences_divisees[i][0];
	}
}

double trouve_pa()