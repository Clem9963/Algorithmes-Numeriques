#include <stdio.h>
#include <stdlib.h>

// PROTOTYPES

double lagrange(double *x, double *y, int n, double a);
double** calcule_differences_divisees(double *x, double *y, int n);
void afficher_differences_divisees(int n, double **differences_divisees);
void liberer_differences_divisees(int n, double **differences_divisees);
double newton(double *x, int n, double **differences_divisees, double a);

int main()
{
	double res = 0;
	/*double x[20] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38};
	double y[20] = {0.99987, 0.99997, 1, 0.99997, 0.99988, 0.99973, 0.99953, 0.99927, 0.99897, 0.99846,
					0.99805, 0.99751, 0.99705, 0.99650, 0.99664, 0.99533, 0.99472, 0.99472, 0.99333, 0.99326};*/

	double x[4] = {1, 2, 3, 4};
	double y[4] = {0, 0, 0, 6};

	double** differences_divisees = calcule_differences_divisees(x, y, 4);
	res = newton(x, 4, differences_divisees, 3);
	printf("%.3f\n", res);
	liberer_differences_divisees(4, differences_divisees);

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

double** calcule_differences_divisees(double *x, double *y, int n)
{
	int i = 0;
	int j = 0;
	double** differences_divisees = calloc(n, sizeof(double*));
	for (i = 0; i < n; i++)
	{
		differences_divisees[i] = calloc(n+1, sizeof(double));
	}

	for (i = 0; i < n; i++)
	{
		differences_divisees[i][0] = x[i];

	}
	for (i = 0; i < n; i++)
	{
		differences_divisees[i][1] = y[i];
	}

	for (j = 2; j < n+1; j++)
	{
		for (i = j-1; i < n; i++)
		{
			differences_divisees[i][j] = (differences_divisees[i][j-1] - differences_divisees[j-2][j-1]) / (differences_divisees[i][0] - differences_divisees[j-2][0]);
		} 
	}
	return differences_divisees;
}

void afficher_differences_divisees(int n, double **differences_divisees)
{
	int i = 0;
	int j = 0;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n+1; j++)
		{
			printf("%.3f   ", differences_divisees[i][j]);
		}
		printf("\n");
	}
}

void liberer_differences_divisees(int n, double **differences_divisees)
{
	int i = 0;
	for (i = 0; i < n; i++)
	{
		free(differences_divisees[i]);
	}
	free(differences_divisees);
}

double newton(double *x, int n, double **differences_divisees, double a)
{
	double* tab_a = calloc(n, sizeof(double));
	double pa = 0;
	double produit_polynomes = 1;
	int i = 0;
	int j = 1;
	
	for (i = 0; i < n; i++)
	{
		tab_a[i] = differences_divisees[i][j];
		j++;
	}

	j = 0;

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < i; j++)
		{
			produit_polynomes *= (a - x[j]);
		}
		pa = pa + (tab_a[i])*produit_polynomes;
		produit_polynomes = 1;
	}

	return pa;
}