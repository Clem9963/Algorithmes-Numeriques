#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void methode_puissance(int ordre, double **A, double *v, double e);
double calcule_norme(double *u, double *v, int ordre);
void multiplier_matrices(double **A, double **B, double **R, int n);
void multiplier_matrice_scalaire(double **A, double k, double **R, int n);

double** allouer_memoire_matrice(int n);
void liberer_memoire_matrice(int n, double **mat);
void afficher_matrice(int n, double **mat);

int main()
{
	double **A = allouer_memoire_matrice(3);
	double **B = allouer_memoire_matrice(3);
	double **R = allouer_memoire_matrice(3);

	A[0][0] = 1;
	A[0][1] = 2;
	A[0][2] = 3;
	A[1][0] = 4;
	A[1][1] = 5;
	A[1][2] = 6;
	A[2][0] = 7;
	A[2][1] = 8;
	A[2][2] = 9;
	
	B[0][0] = 9;
	B[0][1] = 8;
	B[0][2] = 7;
	B[1][0] = 6;
	B[1][1] = 5;
	B[1][2] = 4;
	B[2][0] = 3;
	B[2][1] = 2;
	B[2][2] = 1;

	multiplier_matrices(A, B, R, 3);

	afficher_matrice(3, A);
	afficher_matrice(3, B);
	afficher_matrice(3, R);

	return EXIT_SUCCESS;
}

void methode_puissance(int ordre, double **A, double *v, double e)
{
	// Initialisation de v avec toutes les composantes à 0
	double *v_ancien = calloc(ordre, sizeof(double));
	int i = 0;
	int j = 0;
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

void multiplier_matrices(double **A, double **B, double **R, int n)
{
	int i = 0;
	int j = 0;
	int k = 0;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			R[i][j] = 0;
			for (k = 0; k < n; k++)
			{
				R[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}

void multiplier_matrice_scalaire(double **A, double k, double **R, int n)
{
	int i = 0;
	int j = 0;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			R[i][j] = k * A[i][j];
		}
	}
}

double** allouer_memoire_matrice(int n)
{
	int i = 0;
	double **mat;

	mat = calloc(n, sizeof(double*));
	if (mat == NULL)
	{
		printf("Une erreur est survenue lors de l'allocation de mémoire\n");
		exit(EXIT_FAILURE);
	}

	for (i = 0 ; i < n ; i++)
	{
		mat[i] = calloc(n, sizeof(double));
		if (mat == NULL)
		{
			printf("Une erreur est survenue lors de l'allocation de mémoire\n");
			exit(EXIT_FAILURE);
		}
	}
	return mat;
}

void liberer_memoire_matrice(int n, double **mat)
{
	int i = 0;

	for (i = 0 ; i < n ; i++)
	{
		free(mat[i]);
	}

	free(mat);
}

void afficher_matrice(int n, double **mat)
{
	int i = 0;
	int j = 0;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf("%f\t", mat[i][j]);
		}
		printf("\n");
	}
}
