#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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

/* Fonction principale */
void demarrer_methode(int choix, double *x, double *y, int taille, double borne_inf, double borne_sup);

/* Fonctions d'interpolation */
void lagrange(double *x, double *y, int n, double *a, double *pa, int n_a, opt *mesure);
void newton(double *x, double *y, int n, double *a, double *pa, int n_a, opt *mesure);

/* Fonctions Annexes */
double** calcule_differences_divisees(double *x, double *y, int n, opt *mesure);
void afficher_differences_divisees(int n, double **differences_divisees);
void liberer_differences_divisees(int n, double **differences_divisees);
double* generer_tab_a(double borne_inf, double borne_sup, int *n);
void liberer_tab_a(double *tab_a);


int main()
{
	/*double x[20] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38};
	double y[20] = {0.99987, 0.99997, 1, 0.99997, 0.99988, 0.99973, 0.99953, 0.99927, 0.99897, 0.99846,
					0.99805, 0.99751, 0.99705, 0.99650, 0.99664, 0.99533, 0.99472, 0.99472, 0.99333, 0.99326};*/

	/*double x[11] = {10, 8, 13, 9, 11, 14, 6, 4, 12, 7, 5};
	double y[11] = {9.14, 8.14, 8.74, 8.77, 9.26, 8.1, 6.13, 3.1, 9.13, 7.26, 4.74};*/

	/*double x[11] = {10, 8, 13, 9, 11, 14, 6, 4, 12, 7, 5};
	double y[11] = {7.46, 6.77, 12.74, 7.11, 7.81, 8.84, 6.08, 5.39, 8.15, 6.42, 5.73};*/

	/*double x[21] = {752, 855, 871, 734, 610, 582, 921, 492, 569, 462, 907, 643, 862, 524, 679, 902, 918, 828, 875, 809, 894};
	double y[21] = {85, 83, 162, 79, 81, 83, 281, 81, 81, 80, 243, 84, 84, 82, 226, 260, 82, 186, 77, 223};*/

	/*double x[4] = {1,2,3,4};
	double y[4] = {0,0,0,6};*/

	double x[6] = {-2, 1, 4, -1, 3, -4};
	double y[6] = {-1, 2, 59, 4, 24, -53};

	demarrer_methode(1, x, y, 6, -4.5, 4.5);
    demarrer_methode(2, x, y, 6, -4.5, 4.5);

	return 0;
}

void demarrer_methode(int choix, double *x, double *y, int taille, double borne_inf, double borne_sup)
{
	int i = 0;
	int n = 0;
	FILE *f_dest;

	opt mesure = {0, 0};

	double *a = generer_tab_a(borne_inf, borne_sup, &n);
	
	double* pa = calloc(n, sizeof(double));

	if (choix == 1)
	{
		lagrange(x, y, taille, a, pa, n, &mesure);
		printf("Vous avez choisi d'utiliser la méthode de Lagrange\n");
	}
	else
	{
		newton(x, y, taille, a, pa, n, &mesure);
		printf("Vous avez choisi d'utiliser la méthode de Newton\n");
	}

	f_dest = fopen("polynome.pol", "w");
	for (i = 0; i < n-1; i++)
	{
		fprintf(f_dest, "%f,%f\n", a[i], pa[i]);
	}
	fprintf(f_dest, "%f,%f", a[i], pa[i]);


	printf("La méthode a nécessité %ld clocks processeur et %ld octets ont été alloués\n", mesure.nb_clocks, mesure.octets);

	liberer_tab_a(a);
	free(pa);
}

void lagrange(double *x, double *y, int n, double *a, double *pa, int n_a, opt *mesure)
{
	mesure->nb_clocks = clock();
	int k = 0; mesure->octets += sizeof(int);
	int i = 0; mesure->octets += sizeof(int);
	int j = 0; mesure->octets += sizeof(int);
	double li = 1; mesure->octets += sizeof(double);
	double somme = 0; mesure->octets += sizeof(double);

	for (k = 0; k < n_a; k++)
	{
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				if (j != i)
				{
					li = li * ((a[k]-x[j])/(x[i]-x[j]));
				}
			}
			somme = somme + y[i] * li;
			li = 1;
		}
		pa[k] = somme;
		somme = 0;
	}

	mesure->nb_clocks = clock() - mesure->nb_clocks;
}

void newton(double *x, double *y, int n, double *a, double *pa, int n_a, opt *mesure)
{
	mesure->nb_clocks = clock();
	double** differences_divisees = calcule_differences_divisees(x, y, n, mesure);
	double* tab_a = calloc(n, sizeof(double)); mesure->octets += n*sizeof(double);
	double produit_polynomes = 1; mesure->octets += sizeof(double);
	int i = 0; mesure->octets += sizeof(int);
	int j = 1; mesure->octets += sizeof(int);
	int k = 0; mesure->octets += sizeof(int);

	
	for (i = 0; i < n; i++)
	{
		tab_a[i] = differences_divisees[i][j];
		j++;
	}

	j = 0;

	for (k = 0; k < n_a; k++)
	{
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < i; j++)
			{
				produit_polynomes *= (a[k] - x[j]);
			}
			pa[k] = pa[k] + (tab_a[i])*produit_polynomes;
			produit_polynomes = 1;
		}
	}
	liberer_differences_divisees(n, differences_divisees);
	mesure->nb_clocks = clock() - mesure->nb_clocks;
}

double** calcule_differences_divisees(double *x, double *y, int n, opt *mesure)
{
	int i = 0; mesure->octets += sizeof(int);
	int j = 0; mesure->octets += sizeof(int);
	double** differences_divisees = calloc(n, sizeof(double*)); mesure->octets += n*sizeof(double*);
	for (i = 0; i < n; i++)
	{
		differences_divisees[i] = calloc(n+1, sizeof(double)); mesure->octets += (n+1)*sizeof(double);
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

double* generer_tab_a(double borne_inf, double borne_sup, int *n)
{
	int i = 0;
	double differences = borne_sup - borne_inf;

	*n = (int)(differences / 0.001 + 1);

	double *tab_a = calloc(*n, sizeof(double));

	for (i = 0; i <*n; i++)
	{
		tab_a[i] = borne_inf + i*0.001;
	}

	return tab_a;
}

void liberer_tab_a(double *tab_a)
{
	free(tab_a);
}
