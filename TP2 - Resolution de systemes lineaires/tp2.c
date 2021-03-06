#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define TRUE 1
#define FALSE 0

/***********/
/* TYPEDEF */
/***********/

struct optimisation
{
	clock_t nb_clocks;
	size_t octets;
};
typedef struct optimisation opt;

/**************/
/* PROTOTYPES */
/**************/

/* Fonctions diverses, allocations, affichages */
int min(int a, int b);
double** allouer_memoire_matrice(int n);
void liberer_memoire_matrice(int n, double **mat);
void afficher_matrice_et_second_membre(int n, double **mat, double *second_membre);
void afficher_solutions_chaque_methode(int n, double *x1, double *x2, double *x3, double *x4);

/* Fonctions de résolution de systèmes */
opt gauss(double **A, double *b, double *x, int n);
opt cholesky(double **mat, double *b, double *x, int n);
opt jacobi(double **mat, double *second_membre, double *x, double e, int n, int max_it);
opt gauss_seidel(double **mat, double *second_membre, double *x, double e, int n, int max_it);

/* Fonctions annexes */
int trouver_decomposition(double **mat, double	**r, double **rt, int n);
void resyst_tri_inf(double **mat, double *b1, double *x1, int n);
void resyst_tri_sup(double **mat, double *b1, double *x1, int n);

/* Fonction de test */
void tester_optimisation(double **mat, double *second_membre, double *x, double e, int n, int max_it);
int symetrique(double **mat, int n);
int diagonale_dominante(double **mat, int n);

/* Fonctions de génération de matrices types */
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
	mat = generer_matrice_kms(3, 0.5);

	double second_membre[3];
	double x[3];

    second_membre[0] = 12;
    second_membre[1] = 48;
    second_membre[2]= 57;

    x[0]=0;
    x[1]=0;
    x[2]=0;

	tester_optimisation(mat, second_membre, x,0.001,3,1024);

	liberer_memoire_matrice(3, mat);

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

void afficher_matrice_et_second_membre(int n, double **mat, double *second_membre)
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

void afficher_solutions_chaque_methode(int n, double *x_gauss, double *x_cholesky, double *x_jacobi, double *x_gauss_seidel)
{
	int i = 0;
	printf("Affichage des solutions trouvées par les méthodes suivantes de gauche à droite :\n");
	printf("Gauss, Cholesky, Jacobi, Gauss-Seidel\n");

	for (i = 0 ; i < n ; i++)
	{
		printf("%f\t", x_gauss[i]);
		printf("%f\t", x_cholesky[i]);
		printf("%f\t", x_jacobi[i]);
		printf("%f\t\n", x_gauss_seidel[i]);
	}
}

opt gauss(double **A, double *b, double *x, int n)
{
	opt mesure = {clock(), 0};

	int i, j, k; mesure.octets += 3*sizeof(int);
	int imin; mesure.octets += sizeof(int);
	double p;  mesure.octets += sizeof(double);
	double sum, valmin, tump1, tump2; mesure.octets += 4*sizeof(double);
	
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

	mesure.nb_clocks = clock() - mesure.nb_clocks;
	return mesure;
}

opt cholesky(double **mat, double *b, double *x, int n)
{
	opt mesure = {clock(), 0};

	double *y; mesure.octets += sizeof(double*);
	double **r; mesure.octets += sizeof(double**);
	double **rt; mesure.octets += sizeof(double**);
	r = allouer_memoire_matrice(n); mesure.octets = mesure.octets + n*sizeof(double*) + n*sizeof(double);
	rt = allouer_memoire_matrice(n); mesure.octets = mesure.octets + n*sizeof(double*) + n*sizeof(double);

	y = calloc(n, sizeof(double)); mesure.octets += n*sizeof(double);
	if (y == NULL)
	{
		printf("Une erreur est survenue lors de l'allocation de mémoire\n");
		exit(EXIT_FAILURE);
	}

	if (trouver_decomposition(mat, r, rt, n) ==  FALSE)
	{
		printf("La matrice n'est pas défine positive\n");
		return mesure;
	}

	resyst_tri_inf(rt, b, y, n);
	resyst_tri_sup(r, y, x, n);

	free(y);
	liberer_memoire_matrice(n, r);
	liberer_memoire_matrice(n, rt);

	mesure.nb_clocks = clock() - mesure.nb_clocks;
	return mesure;
}

opt jacobi(double **mat, double *second_membre, double *x, double e, int n, int max_it)
{
	opt mesure = {clock(), 0};

	int i = 0; mesure.octets += sizeof(int);
	int j = 0; mesure.octets += sizeof(int);
	int compteur = 0; mesure.octets += sizeof(int);
	double somme = 0; mesure.octets += sizeof(double);
	double norme = e; mesure.octets += sizeof(double);

	double *y = calloc(n, sizeof(double)); mesure.octets += n*sizeof(double);

	while(norme >= e && compteur < max_it)
	{
		for (i = 0 ; i < n ; i++)
		{
			somme = 0;
			for(j = 0 ; j < i ; j++)
			{
					somme += mat[i][j] * x[j];
			}
			for(j = i + 1 ; j < n ; j++)
			{
					somme += mat[i][j] * x[j];
			}
			y[i] = (second_membre[i] - somme) / mat[i][i];
		}
		somme = 0;

		for (i = 0 ; i < n ; i++){
			somme += pow((y[i] - x[i]), 2);
		}

		for(i = 0 ; i < n ; i++){
			x[i] = y[i];
		}
		norme = sqrtf(somme);
		compteur ++;
	}

	free(y);

	printf("Jacobi : Il y a eu %d itération(s) pour arriver au résultat\n\n", compteur);

	mesure.nb_clocks = clock() - mesure.nb_clocks;
	return mesure;
}

opt gauss_seidel(double **mat, double *second_membre, double *x, double e, int n, int max_it)
{
	opt mesure = {clock(), 0};

	int i = 0; mesure.octets += sizeof(int);
	int j = 0; mesure.octets += sizeof(int);
	int compteur = 0; mesure.octets += sizeof(int);
	double somme = 0; mesure.octets += sizeof(double);
	double norme = e; mesure.octets += sizeof(double);

	double *y = calloc(n, sizeof(double)); mesure.octets += n*sizeof(double);

	while(norme >= e && compteur < max_it)
	{
		for (i = 0 ; i < n ; i++)
		{
			y[i] = x[i];
			somme = 0;
			for(j = 0 ; j < i ; j++)
			{
				somme += mat[i][j] * x[j];
			}
			for(j = i + 1 ; j < n ; j++)
			{
				somme += mat[i][j] * x[j];
			}
			x[i] = (second_membre[i] - somme) / mat[i][i];
		}
		somme = 0;

		for (i = 0 ; i < n ; i++){
			somme += pow((y[i] - x[i]), 2);
		}

		norme = sqrtf(somme);
		compteur ++;
	}

	free(y);

	printf("Gauss-Seidel : Il y a eu %d itération(s) pour arriver au résultat\n\n", compteur);

	mesure.nb_clocks = clock() - mesure.nb_clocks;
	return mesure;
}

int trouver_decomposition(double **mat, double	**r, double **rt, int n)
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
	for (i = 0 ; i < n ; i++)
	{
		for (j = 0 ; j < n ; j++)
		{
			rt[i][j] = r[j][i];
		}
	}
	return TRUE;						// La matrice mat est bien définie positive et décomposable et l'on a trouvé sa décomposition
}

/*****Resolution du systeme lineaire (triangle inferieur)*****/
void resyst_tri_inf(double **mat, double *b1, double *x1, int n)
{
 	int i, k;
	double S;
 	x1[0] = b1[0] / mat[0][0];
 	for (i = 1 ; i < n ; i++)
 	{
		S = 0;
  		for (k = 0 ; k <= i ; k++)
			S += mat[i][k] * x1[k];
		x1[i] = (b1[i]-S) / mat[i][i];
 	}
}

/*****Resolution du systeme lineaire (triangle superieur)*****/
void resyst_tri_sup(double **mat, double *b1, double *x1, int n)
{
 	int i,k;
	double S;
 	x1[n-1] = b1[n-1] / mat[n-1][n-1];
 	for (i = n-2 ; i >= 0 ; i--)
 	{
		S = 0;
  		for (k = n-1 ; k >= i+1 ; k--)
			S += mat[i][k] * x1[k];
		x1[i] = (b1[i]-S) / mat[i][i];
 	}
}

void tester_optimisation(double **mat, double *second_membre, double *x, double e, int n, int max_it)
{
	/* Cette fonction n'a de sens que si la matrice mat convient pour toutes les méthodes : elle doit donc être :
	-> symétrique
	-> définie positive
	-> à diagonale strictement dominante */
	
	if (!symetrique(mat, n) || !diagonale_dominante(mat, n))
	{
		printf("La matrice fournie en argument ne convient pas :\nElle n'est pas symétrique ou non à diagonale dominante\n");
		return;
	}

	int i = 0;
	int j = 0;
	
	opt mesure_gauss = {0, 0};
	opt mesure_cholesky = {0, 0};
	opt mesure_jacobi = {0, 0};
	opt mesure_gauss_seidel = {0, 0};

	double **mat_copy = allouer_memoire_matrice(n);	// Préparation d'une copie de la matrice car elle risque d'être modifiée
	double *x_copy = calloc(n, sizeof(double));
	double *x_gauss = calloc(n, sizeof(double));
	double *x_cholesky = calloc(n, sizeof(double));
	double *x_jacobi = calloc(n, sizeof(double));
	double *x_gauss_seidel = calloc(n, sizeof(double));
	double *second_membre_copy = calloc(n, sizeof(double));

	void reinitialiser()
	{
		for (i = 0 ; i < n ; i++)
		{
			x[i] = x_copy[i];
		}
		for (i = 0 ; i < n ; i++)
		{
			second_membre[i] = second_membre_copy[i];
		}
		for (i = 0 ; i < n ; i++)
		{
			for (j = 0 ; j < n ; j++)
			{
				mat[i][j] = mat_copy[i][j];
			}
		}
	}

	for (i = 0 ; i < n ; i++)
	{
		x_copy[i] = x[i];
	}
	for (i = 0 ; i < n ; i++)
	{
		second_membre_copy[i] = second_membre[i];
	}
	for (i = 0 ; i < n ; i++)
	{
		for (j = 0 ; j < n ; j++)
		{
			mat_copy[i][j] = mat[i][j];
		}
	}

	mesure_gauss = gauss(mat, second_membre, x, n);
	for (i = 0 ; i < n ; i++)
	{
		x_gauss[i] = x[i];
	}
	reinitialiser();

	mesure_cholesky = cholesky(mat, second_membre, x, n);
	for (i = 0 ; i < n ; i++)
	{
		x_cholesky[i] = x[i];
	}
	reinitialiser();

	mesure_jacobi = jacobi(mat, second_membre, x, e, n, max_it);
	for (i = 0 ; i < n ; i++)
	{
		x_jacobi[i] = x[i];
	}
	reinitialiser();

	mesure_gauss_seidel = gauss_seidel(mat, second_membre, x, e, n, max_it);
	for (i = 0 ; i < n ; i++)
	{
		x_gauss_seidel[i] = x[i];
	}

	printf("\tDonnées brutes :\n\n");
	printf("Gauss : %ld clocks processeur, %ld octets alloués\n", mesure_gauss.nb_clocks, mesure_gauss.octets);
	printf("Cholesky : %ld clocks processeur, %ld octets alloués\n", mesure_cholesky.nb_clocks, mesure_cholesky.octets);
	printf("Jacobi : %ld clocks processeur, %ld octets alloués\n", mesure_jacobi.nb_clocks, mesure_jacobi.octets);
	printf("Gauss-Seidel : %ld clocks processeur, %ld octets alloués\n\n\n", mesure_gauss_seidel.nb_clocks, mesure_gauss_seidel.octets);

	printf("\tDonnées de temps d'éxécution en pourcentages par rapport à Gauss:\n\n");
	printf("Cholesky a nécessité %.1f%% de clocks\n", (100*mesure_cholesky.nb_clocks/(double)mesure_gauss.nb_clocks));
	printf("jacobi a nécessité %.1f%% de clocks\n", (100*mesure_jacobi.nb_clocks/(double)mesure_gauss.nb_clocks));
	printf("Gauss-Seidel a nécessité %.1f%% de clocks\n\n\n", (100*mesure_gauss_seidel.nb_clocks/(double)mesure_gauss.nb_clocks));

	printf("\tDonnées d'allocation de mémoire en pourcentages par rapport à Gauss:\n\n");
	printf("Cholesky a nécessité %.1f%% d'octets\n", (100*mesure_cholesky.octets/(double)mesure_gauss.octets));
	printf("jacobi a nécessité %.1f%% d'octets\n", (100*mesure_jacobi.octets/(double)mesure_gauss.octets));
	printf("Gauss-Seidel a nécessité %.1f%% d'octets\n\n\n", (100*mesure_gauss_seidel.octets/(double)mesure_gauss.octets));

	afficher_solutions_chaque_methode(n, x_gauss, x_cholesky, x_jacobi, x_gauss_seidel);
}

int symetrique(double **mat,int n)
{
    int i = 0;
    int j = 0;
    for(i=0 ; i<n ; i++)
    {
        for(j=0 ; j<n ; j++)
        {
            if(mat[i][j] != mat[j][i])
            {
                return FALSE;
            }
        }
    }
    return TRUE;
}

int diagonale_dominante(double **mat,int n)
{
    int i = 0;
    int j = 0;
    int sommeligne = 0;
    int valeur_diago = 0;
    for(i=0 ; i<n ; i++)
    {
        for(j=0 ; j<n ; j++)
        {
            if(i == j)
            {
                valeur_diago = mat[i][j];
            }
            sommeligne += mat[i][j];
        }
        if(fabs(valeur_diago) <= fabs(sommeligne - valeur_diago))
        {
            return FALSE; 
        }
        valeur_diago = 0;
        sommeligne = 0;
    }
    return TRUE;
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
