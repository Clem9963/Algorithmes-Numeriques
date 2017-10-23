#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define TRUE 1
#define FALSE 0

/* PROTOTYPES */

int min(int a, int b);
double** allouer_memoire_matrice(int n);
void liberer_memoire_matrice(int n, double **mat);
void afficher_matrice(int n, double **mat);

void gauss(int n, double **mat, double *second_membre);
void jacobi(double **mat, double *second_membre, double *x, double e, int n, int max_it);

double** generer_matrice_creuse(int n);
double** generer_matrice_bord(int n);
double** generer_matrice_ding_dong(int n);
double** generer_matrice_de_franc(int n);
double** generer_matrice_de_hilbert(int n);
double** generer_matrice_kms(int n, double p);
double** generer_matrice_de_lotkin(int n);
double** generer_matrice_de_moler(int n);

int main()
{
	double **mat;

	mat = generer_matrice_de_moler(3);

	mat[0][0] = 1;
	mat[0][1] = 2;
	mat[0][2] = 2;
	mat[1][0] = 1;
	mat[1][1] = 3;
	mat[1][2] = 3;
	mat[2][0] = 3;
	mat[2][1] = 7;
	mat[2][2] = 8;

	afficher_matrice(3, mat);
	printf("\n\n");

	double second_membre[3];
	second_membre[0] = 2;
	second_membre[1] = 2;
	second_membre[1] = 8;

	double x[3];
	x[0] = 0;
	x[1] = 0;
	x[2] = 0;

	jacobi(mat, second_membre, x, 0.01, 3, 1024);

	printf("%f\n", x[0]);
	printf("%f\n", x[1]);
	printf("%f\n", x[2]);

	return EXIT_SUCCESS;
}

int min(int a, int b)
{
	if (a > b)
	{
		return b;
	}
	else
	{
		return a;
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
	int x = 0;
	int nb_caracteres = 0;

	for (i = 0 ; i < n ; i++)
	{
		for (j = 0 ; j < n ; j++)
		{
			nb_caracteres = printf("%f", mat[i][j]);
			for (x = 0 ; x < 20 - nb_caracteres ; x++)
			{
				printf(" ");
			}
		}
		printf("\n");
	}
}

void gauss(int n, double **mat, double *second_membre)
{
	int i = 0;
	int j = 0;
	int k = 0;
	double premier_coeff = 0;
	
	/*double alpha = 0;

	for (k = 0 ; k < n-1 ; k++)
	{
		for (i = k + 1 ; i < n ; i++)
		{
			alpha = mat[i][k]/mat[k][k];
			for (j = k + 1 ; j < n ; j++)
			{
				mat[i][j] = mat[i][j] - (alpha * second_membre[k]);
				second_membre[i] = second_membre[i] - (alpha * second_membre[k]);
			}
		}
	}*/

	for (j = 0 ; j < n-1 ; j++)
	{
		for (i = j + 1 ; i < n ; i++)
		{
			second_membre[i] = second_membre[i] - (mat[i][j]/mat[j][j]) * second_membre[j];
			for (k = j + 1 ; k < n ; k++)
			{
				mat[i][k] = mat[i][k] - (mat[i][j]/mat[j][j]) * mat[j][k];
			}
			mat[i][j] = 0;
		}
	}
	
	for (i = 0 ; i < n ; i++)
	{
		premier_coeff = mat[i][i];
		for (j =  i ; j < n ; j++)
			{
				mat[i][j] /= premier_coeff;
			}
		second_membre[i] /= premier_coeff;
	}	
}

void jacobi(double **mat, double *second_membre, double *x, double e, int n, int max_it)
{
	int i = 0;
	int j = 0;
	int compteur = 0;
	double somme_avant_i = 0;
	double somme_apres_i = 0;
	int converge = FALSE;

	double *y = calloc(n, sizeof(double));

	while (converge == FALSE)
	{ 
		for (i = 0 ; i < n ; i++)
		{
			for (j = 0 ; j < i ; j++)
			{
				somme_avant_i = somme_avant_i + mat[i][j]*x[j];
			}
			for (j = i + 1 ; j < n ; j++)
			{
				somme_apres_i = somme_apres_i + mat[i][j]*x[j];
			}
			y[i] = second_membre[i] - somme_avant_i - somme_apres_i;
			if (fabs(x[i] - y[i]) < e)
			{
				x[i] = y[i];
				y[i] = 1; 		// On stocke un 1 pour indiquer que la différence entre x[i] et l'ancienne valeur x[i] est inférieure à e
			}
			else
			{
				x[i] = y[i];
				y[i] = 0; 		// 0 sinon
			}
		}
		converge = TRUE;		// Préparation de converge pour le test
		for (i = 0 ; i < n ; i++)
		{
			converge = converge && y[i];
		}
		compteur++;
		converge = converge || (compteur > max_it);
	}
}

double** generer_matrice_creuse(int n)
{
	double **mat;
	int nb_valeurs = 0;
	int i = 0;
	int j = 0;
	int x = 0;
	double al = 0;

	srand(time(NULL));

	mat = allouer_memoire_matrice(n);
	
	nb_valeurs = (n*n*30)/100;
	
	for (x = 0 ; x < nb_valeurs ; x++)
	{
		i = rand() % n;
		j = rand() % n;

		if (mat[i][j] == 0)
		{
			al = (float)rand()/RAND_MAX;
			al -= 0.5;
			al *= M_PI;
			al = tan(al);
			mat[i][j] = al;
		}
		else
		{
			x--;
		}
	}
	return mat;
}

double** generer_matrice_bord(int n)
{
	double **mat;
	int i = 0;
	int j = 0;
	int x = 0;

	srand(time(NULL));
	x = (rand() % 2);

	mat = allouer_memoire_matrice(n);

	for (i = 0 ; i < n ; i++)
	{
		mat[i][i] = 1;
	}

	if (x == 0)
	{
		for (j = 0 ; j < n ; j++)
		{
			mat[0][j] = pow(2, -j);
			mat[j][0] = mat[0][j];
		}
	}
	else
	{
		for (j = 0 ; j < n ; j++)
		{
			mat[n-1][j] = pow(2, n-j-1);
			mat[j][n-1] = mat[n-1][j];
		}
	}

	return mat;
}

double** generer_matrice_ding_dong(int n)
{
	double **mat;
	int i = 0;
	int j = 0;

	mat = allouer_memoire_matrice(n);

	for (i = 0 ; i < n ; i++)
	{
		for(j = 0 ; j < n ; j++)
		{
			mat[i][j] = 1/(2*(n - i - j - 2 + 1.5));
		}
	}

	return mat;
}

double** generer_matrice_de_franc(int n)
{
	double **mat;
	int i = 0;
	int j = 0;

	mat = allouer_memoire_matrice(n);

	for (i = 0 ; i < n ; i++)
	{
		for(j = 0 ; j < n ; j++)
		{
			if (i >= j+2)
			{
				mat[i][j] = 0;
			}
			else
			{
				mat[i][j] = min(i+1, j+1);
			}
		}
	}

	return mat;
}

double** generer_matrice_de_hilbert(int n)
{
	double **mat;
	int i = 0;
	int j = 0;

	mat = allouer_memoire_matrice(n);

	for (i = 0 ; i < n ; i++)
	{
		for(j = 0 ; j < n ; j++)
		{
			mat[i][j] = 1.0/(i + 1 + j);
		}
	}

	return mat;
}

double** generer_matrice_kms(int n, double p)
{
	if (p <= 0 || p >= 1)
	{
		printf("Le paramètre passé en argument est erroné\n");
		return NULL;
	}

	double **mat;
	int i = 0;
	int j = 0;

	mat = allouer_memoire_matrice(n);

	for (i = 0 ; i < n ; i++)
	{
		for(j = 0 ; j < n ; j++)
		{
			mat[i][j] = pow(p, fabs(i-j));
		}
	}

	return mat;
}

double** generer_matrice_de_lotkin(int n)
{
	double **mat;
	int i = 0;
	int j = 0;

	mat = allouer_memoire_matrice(n);

	for (i = 0 ; i < n ; i++)
	{
		for(j = 0 ; j < n ; j++)
		{
			if (i == 0)
			{
				mat[0][j] = 1;
			}
			else
			{
				mat[i][j] = 1.0/(i + 1 + j);
			}
		}
	}
	return mat;
}


double** generer_matrice_de_moler(int n)
{
	double **mat;
	int i = 0;
	int j = 0;

	mat = allouer_memoire_matrice(n);

	for (i = 0 ; i < n ; i++)
	{
		for(j = 0 ; j < n ; j++)
		{
			if (i == j)
			{
				mat[i][j] = i+1;
			}
			else 
			{
				mat[i][j] = min(i+1 ,j+1) - 2;
			}
		}
	}
	return mat;
}
