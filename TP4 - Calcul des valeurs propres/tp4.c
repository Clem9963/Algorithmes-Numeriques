#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

struct optimisation
{
	clock_t nb_clocks;
	size_t octets;
};
typedef struct optimisation opt;

/* PROTOTYPES */

/* Methodes de calcul des valeurs propres*/
opt methode_leverrier_base(int ordre, double **A, double *coefficients);
opt methode_leverrier_amelioree(int ordre, double **A, double *coefficients);
opt methode_puissance(int ordre, double **A, double *v, double e);

/* Methodes de calcul et de traitement brut des matrices*/
void generer_identite(double **A, int n);
void copier_matrice(double **src, double **dest, int n);
void puissance_matrice(double **A, int p, double **R, int n);
void multiplier_matrices(double **A, double **B, double **R, int n);
void multiplier_matrice_scalaire(double **A, double k, double **R, int n);
void additionner_matrices(double **A, double **B, double **R, int n);
double calcule_norme(double *u, double *v, int ordre);
double calcule_trace(double **A, int n);

/* Methodes de gestion de mémoire et d'affichage des matrices*/
double** allouer_memoire_matrice(int n);
void liberer_memoire_matrice(int n, double **mat);
void afficher_matrice(int n, double **mat);

/* Fonctions de génération de matrices types */
double** generer_matrice_creuse(int n);
double** generer_matrice_bord(int n);
double** generer_matrice_ding_dong(int n);
double** generer_matrice_de_franc(int n);
double** generer_matrice_de_hilbert(int n);
double** generer_matrice_kms(int n, double p);
double** generer_matrice_de_lotkin(int n);
double** generer_matrice_de_moler(int n);

/* Fonction annexe */
int min(int a, int b);

int main()
{
	int ordre = 2;
	double e = pow(10, -5);

	double **A = allouer_memoire_matrice(ordre);

	// A[0][0] = 1;
	// A[0][1] = 0;
	// A[0][2] = 0;
	// A[1][0] = 0;
	// A[1][1] = 2;
	// A[1][2] = 0;
	// A[2][0] = 0;
	// A[2][1] = 0;
	// A[2][2] = 3;

	A[0][0] = 1;
	A[0][1] = 0;
	A[1][0] = 0;
	A[1][1] = 2;

	/**/
	double v[ordre];
	opt mesure = {0, 0};
	mesure = methode_puissance(ordre, A, v, e);
	printf("Méthode des puissances : %ld clocks processeur, %ld octets alloués\n", mesure.nb_clocks, mesure.octets);
	/**/

	/*
	int i = 0;
	double coefficients_base[ordre+1];
	double coefficients_amelioree[ordre+1];
	opt mesure_base = {0, 0};
	opt mesure_amelioree = {0, 0};
	mesure_base = methode_leverrier_base(ordre, A, coefficients_base);
	mesure_amelioree = methode_leverrier_amelioree(ordre, A, coefficients_amelioree);
	printf("Coefficients obtenus :\nPar la méthode de Leverrier de base (à gauche), et améliorée (à droite)\n");
	for (i = 0; i < ordre+1; i++)
	{
		printf("%f\t", coefficients_base[i]);
		printf("%f\n", coefficients_amelioree[i]);
	}
	printf("\n");
	printf("Méthode de Leverrier de base : %ld clocks processeur, %ld octets alloués\n", mesure_base.nb_clocks, mesure_base.octets);
	printf("Méthode de Leverrier améliorée : %ld clocks processeur, %ld octets alloués\n", mesure_amelioree.nb_clocks, mesure_amelioree.octets);
	*/

	liberer_memoire_matrice(ordre, A);

	return EXIT_SUCCESS;
}

opt methode_leverrier_base(int ordre, double **A, double *coefficients)
{
	opt mesure = {clock(), 0};

	int i = 0;	mesure.octets += sizeof(int);
	int j = 0;	mesure.octets += sizeof(int);
	double **M = allouer_memoire_matrice(ordre); mesure.octets = mesure.octets + sizeof(double**) + ordre*sizeof(double*) + ordre*sizeof(double);

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
	/* octets alloués par la fonction puissance_matrice */
	mesure.octets = mesure.octets + sizeof(double**) + ordre*sizeof(double*) + ordre*sizeof(double);
	mesure.octets += 6*sizeof(int);

	/* octets alloués par la fonction calcule_trace */
	mesure.octets += sizeof(int);
	mesure.octets += sizeof(double);
	
	liberer_memoire_matrice(ordre, M);

	mesure.nb_clocks = clock() - mesure.nb_clocks;
	return mesure;
}

opt methode_leverrier_amelioree(int ordre, double **A, double *coefficients)
{
	opt mesure = {clock(), 0};

	double **A_buffer = allouer_memoire_matrice(ordre); mesure.octets = mesure.octets + sizeof(double**) + ordre*sizeof(double*) + ordre*sizeof(double);
	copier_matrice(A, A_buffer, ordre);
	double **B_buffer = allouer_memoire_matrice(ordre); mesure.octets = mesure.octets + sizeof(double**) + ordre*sizeof(double*) + ordre*sizeof(double);
	generer_identite(B_buffer, ordre);
	double **M_buffer = allouer_memoire_matrice(ordre); mesure.octets = mesure.octets + sizeof(double**) + ordre*sizeof(double*) + ordre*sizeof(double);
	int i = 0; mesure.octets += sizeof(int);

	coefficients[0] = pow(-1, ordre);

	for (i = 1; i <= ordre; i++)
	{
		multiplier_matrices(B_buffer, A, A_buffer, ordre);
		coefficients[i] = pow(-1, ordre+1)*(1.0/i)*calcule_trace(A_buffer, ordre);
		generer_identite(M_buffer, ordre);
		multiplier_matrice_scalaire(M_buffer, -pow(-1, ordre+1)*coefficients[i], B_buffer, ordre);
		additionner_matrices(A_buffer, B_buffer, B_buffer, ordre);
	}
	/* octets alloués par la fonction multiplier_matrices */
	mesure.octets += 3*sizeof(int);

	/* octets alloués par la fonction calcule_trace */
	mesure.octets += sizeof(int);
	mesure.octets += sizeof(double);

	/* octets alloués par la fonction generer_identite */
	mesure.octets += 2*sizeof(int);

	/* octets alloués par la fonction multiplier_matrice_scalaire */
	mesure.octets += 2*sizeof(int);

	/* octets alloués par la fonction additionner_matrices */
	mesure.octets += 2*sizeof(int);

	liberer_memoire_matrice(ordre, A_buffer);
	liberer_memoire_matrice(ordre, B_buffer);
	liberer_memoire_matrice(ordre, M_buffer);

	mesure.nb_clocks = clock() - mesure.nb_clocks;
	return mesure;
}

opt methode_puissance(int ordre, double **A, double *v, double e)
{
	opt mesure = {clock(), 0};

	// Initialisation de v avec toutes les composantes à 0
	double *v_ancien = calloc(ordre, sizeof(double)); mesure.octets = mesure.octets + sizeof(double*) + ordre*sizeof(double);
	int i = 0; mesure.octets += sizeof(int);
	int j = 0; mesure.octets += sizeof(int);
	double numerateur = 0; mesure.octets += sizeof(double);
	double denominateur = 0; mesure.octets += sizeof(double);

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

		i = ordre-1;
		while (v[i] == 0)		// Traitement du cas où la dernière composante est nulle
		{
			i--;
		}

		for (i = 0; i < ordre; i++)
		{
			v[i] = v[i]/v[ordre-1];
		}
	}

	for (i = 0; i < ordre; i++)		// Test de la convergence
	{
		if (!finite(v[i]))
		{
			printf("La matrice n'admet pas de valeur propre maximale\n");
			mesure.nb_clocks = clock() - mesure.nb_clocks;
			return mesure;
		}
	}

	printf("Voici le vecteur propre associé à la plus grande valeur propre de la matrice\nfournie en argument :\n\n");
	for (i = 0; i < ordre; i++)
	{
		printf("%f\n", v[i]);
	}
	printf("\n");

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

	printf("Voici la plus grande valeur propre de la matrice fournie en argument :\n%f\n", numerateur/denominateur);

	mesure.nb_clocks = clock() - mesure.nb_clocks;
	return mesure;
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
