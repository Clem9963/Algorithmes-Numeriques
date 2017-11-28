#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void methode_puissance(int ordre, double A[ordre][ordre], double *v, double e);
double calcule_norme(double *u, double *v, int ordre);

int main()
{
	double A[2][2] = {{3, 6}, {1, 4}};
	int ordre = 2;
	double v[2] = {0, 0};
	double e = 0.01;

	methode_puissance(ordre, A, v, e);

	return EXIT_SUCCESS;
}

void methode_puissance(int ordre, double A[ordre][ordre], double *v, double e)
{
	// Initialisation de v avec toutes les composantes à 0
	double *v_ancien = calloc(ordre, sizeof(double));
	int i = 0;
	int j = 0;
	int k = 0;
	double numerateur = 0;
	double denominateur = 0;

	printf("%f\n", A[0][0]);
	printf("%f\n", A[0][1]);
	printf("%f\n", A[1][1]);
	printf("%f\n", A[1][0]);

	for (i = 0; i < ordre; i++)
	{
		v[i] = 1;
	}

	while(calcule_norme(v_ancien, v, ordre) > e)
	{
		for (i = 0; i < ordre; i++)
		{
			v_ancien[i] = v[i];
		}

		for (i = 0; i < ordre; i++)
		{
			v[i] = 0;
			for (j = 0; j < ordre; j++)
			{
				v[i] += A[i][j] * v_ancien[j];
			}
		}

		for (i = 0; i < ordre; i++)
		{
			v[i] = v[i]/v[ordre-1];
		}
	}

	printf("Voici le vecteur propre associé à la plus grande valeur propre de la matrice\nfournie en argument :\n");
	for (i = 0; i < ordre; i++)
	{
		printf("%f\n", v[i]);
	}

	for (i = 0; i < ordre; i++)
	{
		v_ancien[i] = 0;
		for (j = 0; j < ordre; j++)
		{
			v_ancien[i] += A[i][j] * v[j];
		}
	}

	for (i = 0; i < ordre; i++)
	{
		numerateur += v[i] * v_ancien[i];
	}

	for (i = 0; i < ordre; i++)
	{
		denominateur += v[i] * v[i];
	}

	printf("Voici la plus grande valeur propre de la matrice\nfournie en argument : %f\n", numerateur/denominateur);

}

double calcule_norme(double *u, double *v, int ordre)
{
	int i = 0;
	double somme = 0;

	for (i = 0 ; i < ordre ; i++)
	{
		somme += pow((u[i] - v[i]), 2);
	}

	return sqrt(somme);
}