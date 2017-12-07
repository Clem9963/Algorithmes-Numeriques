#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void methode_leverrier_base(int ordre, double **A, double *coefficients);
void methode_leverrier_amelioree(int ordre, double **A, double *coefficients);
void methode_puissance(int ordre, double **A, double *v, double e);

void generer_identite(double **A, int n);
void copier_matrice(double **src, double **dest, int n);
void puissance_matrice(double **A, int p, double **R, int n);
void multiplier_matrices(double **A, double **B, double **R, int n);
void multiplier_matrice_scalaire(double **A, double k, double **R, int n);
void additionner_matrices(double **A, double **B, double **R, int n);
double calcule_norme(double *u, double *v, int ordre);
double calcule_trace(double **A, int n);

double** allouer_memoire_matrice(int n);
void liberer_memoire_matrice(int n, double **mat);
void afficher_matrice(int n, double **mat);

int main()
{
	double **A = allouer_memoire_matrice(3);
	double **B = allouer_memoire_matrice(3);
	//double **R = allouer_memoire_matrice(3);

	A[0][0] = 1;
	A[0][1] = 0;
	A[0][2] = 0;
	A[1][0] = 0;
	A[1][1] = 2;
	A[1][2] = 0;
	A[2][0] = 0;
	A[2][1] = 0;
	A[2][2] = 3;

	B[0][0] = 1;
	B[0][1] = 0;
	B[0][2] = 0;
	B[1][0] = 0;
	B[1][1] = 2;
	B[1][2] = 0;
	B[2][0] = 0;
	B[2][1] = 0;
	B[2][2] = 3;

	double coefficients[4];

	// multiplier_matrices(A, B, R, 3);
	// puissance_matrice(A, 3, R, 3);
	methode_leverrier_amelioree(3, A, coefficients);

	// afficher_matrice(3, A);
	//afficher_matrice(3, R);

	printf("%f\n", coefficients[0]);
	printf("%f\n", coefficients[1]);
	printf("%f\n", coefficients[2]);
	printf("%f\n", coefficients[3]);

	return EXIT_SUCCESS;
}

void methode_leverrier_base(int ordre, double **A, double *coefficients)
{
	int i = 0;
	int j = 0;
	double **M = allouer_memoire_matrice(ordre);

	coefficients[0] = pow(-1, ordre);

	for (i = 1; i <= ordre; i++)
	{
		coefficients[i] = 0;
		for (j = 0; j < i; j++)
		{
			puissance_matrice(A, i-j, M, ordre);
			coefficients[i] += coefficients[j]*calcule_trace(M, ordre);
		}
		coefficients[i] = -coefficients[i]/i;
	}
}

void methode_leverrier_amelioree(int ordre, double **A, double *coefficients)
{
	double **A_buffer = allouer_memoire_matrice(ordre);
	copier_matrice(A, A_buffer, ordre);
	double **B_buffer = allouer_memoire_matrice(ordre);
	generer_identite(B_buffer, ordre);
	double **M_buffer = allouer_memoire_matrice(ordre);
	int i = 0;

	coefficients[0] = pow(-1, ordre);

	for (i = 1; i <= ordre; i++)
	{
		multiplier_matrices(B_buffer, A, A_buffer, ordre);
		coefficients[i] = pow(-1, ordre+1)*(1.0/i)*calcule_trace(A_buffer, ordre);
		generer_identite(M_buffer, ordre);
		multiplier_matrice_scalaire(M_buffer, -pow(-1, ordre+1)*coefficients[i], B_buffer, ordre);
		additionner_matrices(A_buffer, B_buffer, B_buffer, ordre);
	}
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

void generer_identite(double **A, int n)
{
	int i = 0;
	int j = 0;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (i == j)
			{
				A[i][j] = 1;
			}
			else
			{
				A[i][j] = 0;
			}
		}
	}
}

void copier_matrice(double **src, double **dest, int n)
{
	int i = 0;
	int j = 0;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			dest[i][j] = src[i][j];
		}
	}
}

void puissance_matrice(double **A, int p, double **R, int n)
{
	double **M = allouer_memoire_matrice(n);
	int i = 0;
	int j = 0;
	int k = 0;

	for (j = 0; j < n; j++)
	{
		for (k = 0; k < n; k++)
		{
			R[j][k] = A[j][k];
		}
	}
	for (i = 2; i <= p; i++)
	{
		for (j = 0; j < n; j++)
		{
			for (k = 0; k < n; k++)
			{
				M[j][k] = R[j][k];
			}
		}
		multiplier_matrices(A, M, R, n);
	}
	liberer_memoire_matrice(n, M);
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

void additionner_matrices(double **A, double **B, double **R, int n)
{
	int i = 0;
	int j = 0;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			R[i][j] = A[i][j] + B[i][j];
		}
	}
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

double calcule_trace(double **A, int n)
{
	int i = 0;
	double trace = 0;
	for (i = 0; i < n; i++)
	{
		trace += A[i][i];
	}
	return trace;
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
