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
void afficher_matrice(int n, double **mat, double *second_membre);
void afficher_solutions(int n, double *x);

void gauss(double **A, double *b, double *x, int n);
void jacobi(double **mat, double *second_membre, double *x, double e, int n, int max_it);
int trouver_decomposition(double **mat, double	**r, int n);
void resyst_tri_inf(double **mat, double *b1, double *x1, int n);
void resyst_tri_sup(double **mat, double *b1, double *x1, int n);

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
	double **r;

	mat = allouer_memoire_matrice(3);
	r = allouer_memoire_matrice(3);

	mat[0][0] = 1;
	mat[0][1] = 1;
	mat[0][2] = 1;
	mat[1][0] = 1;
	mat[1][1] = 2;
	mat[1][2] = 2;
	mat[2][0] = 1;
	mat[2][1] = 2;
	mat[2][2] = 3;

	double second_membre[3];
	second_membre[0] = 1;
	second_membre[1] = 1;
	second_membre[2] = 1;

	afficher_matrice(3, mat, second_membre);
	printf("\n\n");

	trouver_decomposition(mat, r, 3);
	afficher_matrice(3, r, second_membre);
	printf("\n\n");

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

void afficher_matrice(int n, double **mat, double *second_membre)
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
		nb_caracteres = printf("|");
		for (x = 0 ; x < 10 - nb_caracteres ; x++)
		{
			printf(" ");
		}
		nb_caracteres = printf("%f", second_membre[i]);
		printf("\n");
	}
}

void gauss(double **A, double *b, double *x, int n)
{
	int i, j, k;
	int imin;
	double p;
	double sum, valmin, tump1, tump2;
	
	for (k = 0 ; k < n-1 ; k++)
	{
		/* Dans un premier temps, on cherche l'élément minimum (non */
		/* nul) en valeur absolue dans la colonne k et d'indice i	*/
		/* supérieur ou égal à k.									*/
		
		valmin = A[k][k]; imin = k;
		for (i = k+1 ; i < n ; i++)
		{
			if (valmin != 0)
			{
				if (abs(A[i][k]) < abs(valmin) && A[i][k] != 0)
				{
					valmin = A[i][k];
					imin = i;
				}
			}
			else 
			{
				valmin = A[i][k];
				imin = i;
			}	 
		}
		
		/* Si l'élément minimum est nul, on peut en déduire */
		/* que la matrice est singulière. Le pogramme est	*/
		/* alors interrompu.								*/
		
		if (valmin == 0.)
		{
			printf("\n\n\nAttention! Matrice singuliere!\n\n\n");
			exit(EXIT_FAILURE);
		}
		
		/* Si la matrice n'est pas singulière, on inverse	*/
		/* les éléments de la ligne imax avec les éléments	*/
		/* de la ligne k. On fait de même avec le vecteur b. */
		
		for(j = 0 ; j < n ; j++)
		{
			tump1 = A[imin][j];
			A[imin][j] = A[k][j];
			A[k][j] = tump1;
		}
		
		tump2 = b[imin];
		b[imin] = b[k];
		b[k] = tump2;
		
		
		/* On procède à la réduction de la matrice par la */
		/* méthode d'élimination de Gauss. */
		
		for (i = k+1 ; i < n ; i++)
		{
			p = A[i][k]/A[k][k];
			
			for (j = 0 ; j < n ; j++)
			{
				A[i][j] = A[i][j] - p*A[k][j];
			}
			
			b[i] = b[i] - p*b[k]; 
		}
	}	
	 
	/* On vérifie que la matrice n'est toujours pas singulière. */
	/* Si c'est le cas, on interrompt le programme. */
	 
	if (A[n-1][n-1] == 0)
	{
		printf("\n\n\nAttention! Matrice singuliere!\n\n\n");
		exit(EXIT_FAILURE); 
	}
	 
	/* Une fois le système réduit, on obtient une matrice triangulaire */
	/* supérieure et la résolution du système se fait très simplement. */
	 
	x[n-1] = b[n-1]/A[n-1][n-1];
	 
	for(i = n-2 ; i > -1 ; i--)
	{
			sum = 0;
			
			for(j = n-1 ; j > i ; j--)
			{
				sum = sum + A[i][j]*x[j];
			}
			x[i] = (b[i] - sum)/A[i][i];
	}
}

void jacobi(double **mat, double *second_membre, double *x, double e, int n, int max_it)
{
	int i = 0;
	int j = 0;
	int compteur = 0;
	double somme_avant_i = 0;
	double somme_apres_i = 0;
	double norme = 0;
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
			norme = norme + pow((x[i] - y[i]), 2);
			somme_avant_i = 0;
			somme_apres_i = 0;
		}
		norme = sqrt(norme);
		if (norme < e)
		{
			converge = TRUE;
		}
		norme = 0;
		for (i = 0 ; i < n ; i++)
		{
			x[i] = y[i];
		}
		compteur++;
		converge = converge || (compteur > max_it);
	}
}

int trouver_decomposition(double **mat, double	**r, int n)
{
	int i = 0;
	int j = 0;
	int k = 0;
	double s = 0;
	double somme = 0;
	for(i = 0 ; i < n ; i++)
	{
		for (k = 0 ; k < i ; k++)
		{
			somme = somme + pow(r[k][i], 2); 
		}
		s = mat[i][i] - somme;
		if (s <= 0)
		{
			return FALSE;				// On retourne 0 car A n'est pas définie positive */
		}
		else
		{
			r[i][i] = sqrt(s);
			somme = 0;
			for (j = i + 1 ; j < n ; j++)
			{
				for (k = 0 ; k < i ; k++)
				{
					somme = somme + r[k][i] * r[k][j];
				}
				r[i][j] = (mat[i][j] - somme)/r[i][i];
				somme = 0;
			}
		}
	}
	return TRUE;						// La matrice mat est bien définie positive et décomposable et l'on a trouvé sa décomposition
}

/*****Resolution du systeme lineaire (triangle inferieur)*****/
void resyst_tri_inf(double **mat, double *b1, double *x1, int n)
{
 	int i,k;
	double S;
 	x1[1] = b1[1] / mat[1][1];
 	for (i = 2 ; i <= n ; i++)
 	{
		S = 0;
  		for (k = 1 ; k <= i ; k++)
			S += mat[i][k] * x1[k];
		x1[i] = (b1[i]-S) / mat[i][i];
 	}
}

/*****Resolution du systeme lineaire (triangle superieur)*****/
void resyst_tri_sup(double **mat, double *b1, double *x1, int n)
{
 	int i,k;
	double S;
 	x1[n] = b1[n] / mat[n][n];
 	for (i = n-1 ; i>=1 ; i--)
 	{
		S = 0;
  		for (k = n ; k >= i+2 ; k--)
			S += mat[i][k] * x1[k];
		x1[i] = (b1[i]-S) / mat[i][i];
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
