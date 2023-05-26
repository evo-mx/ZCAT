/*
 * zcat-benchmark.c
 *
 *  Created on: Apr 1, 2022
 *      Author: Saul Zapotecas
 */

//#include "../zcat/sz.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "zcat-F.h"
#include "zcat-g.h"
#include "zcat-Z.h"
#include "zcat-benchmark.h"
#include "zcat-tools.h"

#define ZCAT_TRUE 1
#define ZCAT_FALSE 0

/** Name of the test instances ***************************************************************** **/
char ZCAT_MOP_STR[20][7] =
{ "ZCAT1", "ZCAT2", "ZCAT3", "ZCAT4", "ZCAT5", "ZCAT6", "ZCAT7", "ZCAT8",
		"ZCAT9", "ZCAT10", "ZCAT11", "ZCAT12", "ZCAT13", "ZCAT14", "ZCAT15",
		"ZCAT16", "ZCAT17", "ZCAT18", "ZCAT19", "ZCAT20" };

/** Pointer to the test instances ************************************************************** **/
void
(*ZCAT_MOP[20])(double *X, double *F, double *G,
		int *Xbin) =
		{	ZCAT1, ZCAT2, ZCAT3, ZCAT4, ZCAT5, ZCAT6, ZCAT7, ZCAT8, ZCAT9, ZCAT10, ZCAT11, ZCAT12,
			ZCAT13, ZCAT14, ZCAT15, ZCAT16, ZCAT17, ZCAT18, ZCAT19, ZCAT20
};

/** Pointer to the g functions ***************************************************************** **/
double *
(*ZCAT_G[20])(double *y, int m, int n) =
{
	//  F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14, F15, F16, F17, F18, F19, F20
		g4, g5, g2, g7, g9, g4, g5, g2, g7, g9, g3, g10, g1, g6, g8, g10, g1, g8, g6, g3
};

/** Pointer to the F functions ***************************************************************** **/
void
(*ZCAT_F[20])(double *F, double *y,
		int M) =
		{	F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14, F15, F16, F17, F18, F19, F20
};

struct
{
	int M; /* Number of objectives */
	int N; /* Number of decision variables */
	int COMPLICATED_PS; /* Complicated PS flag */
	int LEVEL; /* Difficulty level */
	int BIAS; /* Bias flag */
	int IMBALANCE; /* Imbalance MOP */

	double *LB; /* Low bound */
	double *UB; /* Up bound */
} zcat_config;

/**
 * Set the bounds for all test problems
 * @param LB:   low bound
 * @param UB:   up bound
 * @param n:    number of decision variables
 */
void zcat_set_bounds(double *LB, double *UB, int n)
{
	int i;
	for (i = 1; i < n + 1; ++i)
	{
		LB[i - 1] = -i * 0.5;
		UB[i - 1] = i * 0.5;
	}
	return;
}

/**
 * get the normalization of x
 * @param x:    Decision vector
 * @param LB:   Low bound
 * @param UB:   Up bound
 * @param n:    Number of decision variables
 * @return:     Normalization of x into the range [0,1]
 */
double *
zcat_get_y(double *x, double *LB, double *UB, int n)
{
	int i;
	double *y;

	y = (double*) malloc(sizeof(double) * n);
	for (i = 0; i < n; ++i)
	{
		y[i] = (x[i] - LB[i]) / (UB[i] - LB[i]);
		assert(0.0 <= y[i] && y[i] <= 1.0);
	}
	return (y);
}

/**
 * Define z_{1:n-m} = (y_{m+1}-g_{m+1}, ..., y_n-g_n)
 * @param y:    Decision variables
 * @param m:    Dimension of PF
 * @param n:    Number of decision variables
 * @param g:    g function (g0, g1, ..., g10)
 * @return
 */
double *
zcat_get_z(double *y, int m, int n, double *
(*g_function)(double *y, int m, int n))
{
	int j, i;
	double *z;
	double *g;
	double diff;

	g = g_function(y, m, n);
	z = (double*) malloc(sizeof(double) * (n - m));
	for (i = 0, j = m + 1; j <= n; ++j, ++i)
	{
		assert(i < n - m);
		diff = y[j - 1] - g[i];
		z[i] = (fabs(diff) < DBL_EPSILON) ? 0.0 : diff; /*precision issues*/
		assert(-1.0 <= z[i] && z[i] <= 1.0);
	}

	free(g);
	return z;
}

/**
 * Define w_{1:n-m} = (Zbias(z_{1}),..., Zbias(z_{n-m})) if bias flag is activated,
 * @param z:    z_{1:n-m} = (y_{m+1}-g_{m+1}, ..., y_n-g_n)
 * @param m:    Dimension of PF
 * @param n:    Number of decision variables
 * @return
 */
double *
zcat_get_w(double *z, int m, int n)
{
	int i;
	double *w;

	w = (double*) malloc(sizeof(double) * (n - m));

	for (i = 0; i < n - m; ++i)
	{
		w[i] = (zcat_config.BIAS == 1) ? Zbias(z[i]) : z[i];
	}
	return w;
}

/**
 * Modular approach.
 * @param i:     Objective function
 * @param M:     Number of objectives
 * @param w:     The vector w=YII-g(YI)
 * @param wsize: The size of w (m-n)
 * @param Jsize: The size of J
 * @return
 */
double *
zcat_get_J(int i, int M, double *w, int wsize, int *Jsize)
{
	double *J = NULL;
	int j, size;

	assert(i <= i && i <= M);
	size = 0;
	for (j = 1; j <= wsize; ++j)
	{
		if ((j - i) % M == 0)
		{
			size++;
			J = (double*) realloc(J, sizeof(double) * size);
			J[size - 1] = w[j - 1];
		}
	}
	if (size == 0)
	{
		size = 1;
		J = (double*) realloc(J, sizeof(double) * size);
		J[size - 1] = w[0];
	}
	*Jsize = size;
	assert(size > 0);
	return J;
}

/**
 * Evaluate function Z according to the difficulty level
 * @param w:        The 'w' vector given by the modular approach
 * @param wsize:    The size of 'w'
 * @return
 */
double zcat_evaluate_Z(double *w, int wsize, int ith_objective)
{
	double z;

	if (zcat_config.IMBALANCE == ZCAT_TRUE)
	{
		if (ith_objective % 2 == 0)
		{
			z = Z4(w, wsize);
		}
		else
		{
			z = Z1(w, wsize);
		}
		return z;
	}

	switch (zcat_config.LEVEL)
	{
	/* Basis functions */
	case 1:
		z = Z1(w, wsize);
		break;
	case 2:
		z = Z2(w, wsize);
		break;
	case 3:
		z = Z3(w, wsize);
		break;
	case 4:
		z = Z4(w, wsize);
		break;
	case 5:
		z = Z5(w, wsize);
		break;
	case 6:
		z = Z6(w, wsize);
		break;
	default:
		z = Z1(w, wsize);
		break;
	}
	return z;
}

/**
 * Obtain the values of B's functions
 * @param y:        The y vector (normalized x)
 * @param M:        The number of objectives
 * @param m:        The dimension of PS/PF
 * @param n:        The number of decision variables
 * @return
 */
double *
zcat_get_beta(double *y, int M, int m, int n, double *
(*g_function)(double *y, int m, int n))
{
	double *J, *b, Zvalue;
	int Jsize, i;
	double *w, *z;

	z = zcat_get_z(y, m, n, g_function);
	w = zcat_get_w(z, m, n);

	b = (double*) malloc(sizeof(double) * M);
	// TODO: The following 'if sentnce' was added
	if (n == m) // number of variables equal to the dimension of the PF
	{
		for (i = 1; i <= M; ++i)
		{
			b[i - 1] = 0.0;
		}
	}
	else
	{
		for (i = 1; i <= M; ++i)
		{
			/* Modular approach */
			J = zcat_get_J(i, M, w, n - m, &Jsize);
			Zvalue = zcat_evaluate_Z(J, Jsize, i);
			b[i - 1] = (i * i) * Zvalue;
			free(J);
		}
	}
	free(z);
	free(w);
	return b;
}

/**
 * Obtaining the values of alpha (eq. 5 in the paper)
 * @param y
 * @param M
 * @param f_function
 * @return
 */
double *
zcat_get_alpha(double *y, int M,
		void (*f_function)(double *F, double *y, int M))
{
	int i;
	double *a;
	a = (double*) malloc(sizeof(double) * M);

	f_function(a, y, M);
	for (i = 1; i <= M; ++i)
	{
		a[i - 1] = (i * i) * a[i - 1];
	}
	return a;
}

/**
 * Assigning objective values (additive approach)
 * @param f:        Objectives vector
 * @param alpha:    Alpha values
 * @param beta:     Beta values
 * @param M:        Number of objectives
 */
void zcat_mop_definition(double *f, double *alpha, double *beta, int M)
{
	int i;
	for (i = 1; i <= M; ++i)
	{
		f[i - 1] = alpha[i - 1] + beta[i - 1]; /* Additive Approach */
	}
	return;
}

/** ********************************************************************************************* **
 *  Definition of the problems ZCAT1--ZCAT20
 ** ********************************************************************************************* **/

/**
 * ZCAT1 test function
 * @param x:    Variables vector
 * @param f:    Objectives vector
 * @param g:    Constraints vector (no constraints, input: NULL)
 * @param xbin: Binary vector (no binary vector, input: NULL)
 */
void ZCAT1(double *x, double *f, double *g, int *xbin)
{
	int M = zcat_config.M;
	int n = zcat_config.N;
	int m = zcat_config.M - 1;
	double *y, *alpha, *beta;
	double *
	(*g_function)(double *y, int m, int n);
	void
	(*f_function)(double *F, double *y, int M);

	f_function = ZCAT_F[zcat1];
	g_function = (zcat_config.COMPLICATED_PS) ? ZCAT_G[zcat1] : g0;

	y = zcat_get_y(x, zcat_config.LB, zcat_config.UB, n); /* Normalization */
	alpha = zcat_get_alpha(y, M, f_function); /* Define Position */
	beta = zcat_get_beta(y, M, m, n, g_function); /* Define Distance*/
	zcat_mop_definition(f, alpha, beta, M); /* Assign fitness values */

	free(y);
	free(alpha);
	free(beta);
	return;
}

/**
 * ZCAT2 test function
 * @param x:    Variables vector
 * @param f:    Objectives vector
 * @param g:    Constraints vector (no constraints, input: NULL)
 * @param xbin: Binary vector (no binary vector, input: NULL)
 */
void ZCAT2(double *x, double *f, double *g, int *xbin)
{
	int M = zcat_config.M;
	int n = zcat_config.N;
	int m = zcat_config.M - 1;
	double *y, *alpha, *beta;
	double *
	(*g_function)(double *y, int m, int n);
	void
	(*f_function)(double *F, double *y, int M);

	f_function = ZCAT_F[zcat2];
	g_function = (zcat_config.COMPLICATED_PS) ? ZCAT_G[zcat2] : g0;

	y = zcat_get_y(x, zcat_config.LB, zcat_config.UB, n); /* Normalization */
	alpha = zcat_get_alpha(y, M, f_function); /* Define Position */
	beta = zcat_get_beta(y, M, m, n, g_function); /* Define Distance*/
	zcat_mop_definition(f, alpha, beta, M); /* Assign fitness values */

	free(y);
	free(alpha);
	free(beta);
	return;
}

/**
 * ZCAT3 test function
 * @param x:    Variables vector
 * @param f:    Objectives vector
 * @param g:    Constraints vector (no constraints, input: NULL)
 * @param xbin: Binary vector (no binary vector, input: NULL)
 */
void ZCAT3(double *x, double *f, double *g, int *xbin)
{
	int M = zcat_config.M;
	int n = zcat_config.N;
	int m = zcat_config.M - 1;
	double *y, *alpha, *beta;
	double *
	(*g_function)(double *y, int m, int n);
	void
	(*f_function)(double *F, double *y, int M);

	f_function = ZCAT_F[zcat3];
	g_function = (zcat_config.COMPLICATED_PS) ? ZCAT_G[zcat3] : g0;

	y = zcat_get_y(x, zcat_config.LB, zcat_config.UB, n); /* Normalization */
	alpha = zcat_get_alpha(y, M, f_function); /* Define Position */
	beta = zcat_get_beta(y, M, m, n, g_function); /* Define Distance*/
	zcat_mop_definition(f, alpha, beta, M); /* Assign fitness values */

	free(y);
	free(alpha);
	free(beta);
	return;
}

/**
 * ZCAT4 test function
 * @param x:    Variables vector
 * @param f:    Objectives vector
 * @param g:    Constraints vector (no constraints, input: NULL)
 * @param xbin: Binary vector (no binary vector, input: NULL)
 */
void ZCAT4(double *x, double *f, double *g, int *xbin)
{
	int M = zcat_config.M;
	int n = zcat_config.N;
	int m = zcat_config.M - 1;
	double *y, *alpha, *beta;
	double *
	(*g_function)(double *y, int m, int n);
	void
	(*f_function)(double *F, double *y, int M);

	f_function = ZCAT_F[zcat4];
	g_function = (zcat_config.COMPLICATED_PS) ? ZCAT_G[zcat4] : g0;

	y = zcat_get_y(x, zcat_config.LB, zcat_config.UB, n); /* Normalization */
	alpha = zcat_get_alpha(y, M, f_function); /* Define Position */
	beta = zcat_get_beta(y, M, m, n, g_function); /* Define Distance*/
	zcat_mop_definition(f, alpha, beta, M); /* Assign fitness values */

	free(y);
	free(alpha);
	free(beta);
	return;
}

/**
 * ZCAT5 test function
 * @param x:    Variables vector
 * @param f:    Objectives vector
 * @param g:    Constraints vector (no constraints, input: NULL)
 * @param xbin: Binary vector (no binary vector, input: NULL)
 */
void ZCAT5(double *x, double *f, double *g, int *xbin)
{
	int M = zcat_config.M;
	int n = zcat_config.N;
	int m = zcat_config.M - 1;
	double *y, *alpha, *beta;
	double *
	(*g_function)(double *y, int m, int n);
	void
	(*f_function)(double *F, double *y, int M);

	f_function = ZCAT_F[zcat5];
	g_function = (zcat_config.COMPLICATED_PS) ? ZCAT_G[zcat5] : g0;

	y = zcat_get_y(x, zcat_config.LB, zcat_config.UB, n); /* Normalization */
	alpha = zcat_get_alpha(y, M, f_function); /* Define Position */
	beta = zcat_get_beta(y, M, m, n, g_function); /* Define Distance*/
	zcat_mop_definition(f, alpha, beta, M); /* Assign fitness values */

	free(y);
	free(alpha);
	free(beta);
	return;
}

/**
 * ZCAT6 test function
 * @param x:    Variables vector
 * @param f:    Objectives vector
 * @param g:    Constraints vector (no constraints, input: NULL)
 * @param xbin: Binary vector (no binary vector, input: NULL)
 */
void ZCAT6(double *x, double *f, double *g, int *xbin)
{
	int M = zcat_config.M;
	int n = zcat_config.N;
	int m = zcat_config.M - 1;
	double *y, *alpha, *beta;
	double *
	(*g_function)(double *y, int m, int n);
	void
	(*f_function)(double *F, double *y, int M);

	f_function = ZCAT_F[zcat6];
	g_function = (zcat_config.COMPLICATED_PS) ? ZCAT_G[zcat6] : g0;

	y = zcat_get_y(x, zcat_config.LB, zcat_config.UB, n); /* Normalization */
	alpha = zcat_get_alpha(y, M, f_function); /* Define Position */
	beta = zcat_get_beta(y, M, m, n, g_function); /* Define Distance*/
	zcat_mop_definition(f, alpha, beta, M); /* Assign fitness values */

	free(y);
	free(alpha);
	free(beta);
	return;
}

/**
 * ZCAT7 test function
 * @param x:    Variables vector
 * @param f:    Objectives vector
 * @param g:    Constraints vector (no constraints, input: NULL)
 * @param xbin: Binary vector (no binary vector, input: NULL)
 */
void ZCAT7(double *x, double *f, double *g, int *xbin)
{
	int M = zcat_config.M;
	int n = zcat_config.N;
	int m = zcat_config.M - 1;
	double *y, *alpha, *beta;
	double *
	(*g_function)(double *y, int m, int n);
	void
	(*f_function)(double *F, double *y, int M);

	f_function = ZCAT_F[zcat7];
	g_function = (zcat_config.COMPLICATED_PS) ? ZCAT_G[zcat7] : g0;

	y = zcat_get_y(x, zcat_config.LB, zcat_config.UB, n); /* Normalization */
	alpha = zcat_get_alpha(y, M, f_function); /* Define Position */
	beta = zcat_get_beta(y, M, m, n, g_function); /* Define Distance*/
	zcat_mop_definition(f, alpha, beta, M); /* Assign fitness values */

	free(y);
	free(alpha);
	free(beta);
}

/**
 * ZCAT8 test function
 * @param x:    Variables vector
 * @param f:    Objectives vector
 * @param g:    Constraints vector (no constraints, input: NULL)
 * @param xbin: Binary vector (no binary vector, input: NULL)
 */
void ZCAT8(double *x, double *f, double *g, int *xbin)
{
	int M = zcat_config.M;
	int n = zcat_config.N;
	int m = zcat_config.M - 1;
	double *y, *alpha, *beta;
	double *
	(*g_function)(double *y, int m, int n);
	void
	(*f_function)(double *F, double *y, int M);

	f_function = ZCAT_F[zcat8];
	g_function = (zcat_config.COMPLICATED_PS) ? ZCAT_G[zcat8] : g0;

	y = zcat_get_y(x, zcat_config.LB, zcat_config.UB, n); /* Normalization */
	alpha = zcat_get_alpha(y, M, f_function); /* Define Position */
	beta = zcat_get_beta(y, M, m, n, g_function); /* Define Distance*/
	zcat_mop_definition(f, alpha, beta, M); /* Assign fitness values */

	free(y);
	free(alpha);
	free(beta);
	return;
}

/**
 * ZCAT9 test function
 * @param x:    Variables vector
 * @param f:    Objectives vector
 * @param g:    Constraints vector (no constraints, input: NULL)
 * @param xbin: Binary vector (no binary vector, input: NULL)
 */
void ZCAT9(double *x, double *f, double *g, int *xbin)
{
	int M = zcat_config.M;
	int n = zcat_config.N;
	int m = zcat_config.M - 1;
	double *y, *alpha, *beta;
	double *
	(*g_function)(double *y, int m, int n);
	void
	(*f_function)(double *F, double *y, int M);

	f_function = ZCAT_F[zcat9];
	g_function = (zcat_config.COMPLICATED_PS) ? ZCAT_G[zcat9] : g0;

	y = zcat_get_y(x, zcat_config.LB, zcat_config.UB, n); /* Normalization */
	alpha = zcat_get_alpha(y, M, f_function); /* Define Position */
	beta = zcat_get_beta(y, M, m, n, g_function); /* Define Distance*/
	zcat_mop_definition(f, alpha, beta, M); /* Assign fitness values */

	free(y);
	free(alpha);
	free(beta);
	return;
}

/**
 * ZCAT10 test function
 * @param x:    Variables vector
 * @param f:    Objectives vector
 * @param g:    Constraints vector (no constraints, input: NULL)
 * @param xbin: Binary vector (no binary vector, input: NULL)
 */
void ZCAT10(double *x, double *f, double *g, int *xbin)
{
	int M = zcat_config.M;
	int n = zcat_config.N;
	int m = zcat_config.M - 1;
	double *y, *alpha, *beta;
	double *
	(*g_function)(double *y, int m, int n);
	void
	(*f_function)(double *F, double *y, int M);

	f_function = ZCAT_F[zcat10];
	g_function = (zcat_config.COMPLICATED_PS) ? ZCAT_G[zcat10] : g0;

	y = zcat_get_y(x, zcat_config.LB, zcat_config.UB, n); /* Normalization */
	alpha = zcat_get_alpha(y, M, f_function); /* Define Position */
	beta = zcat_get_beta(y, M, m, n, g_function); /* Define Distance*/
	zcat_mop_definition(f, alpha, beta, M); /* Assign fitness values */

	free(y);
	free(alpha);
	free(beta);
	return;
}

/**
 * ZCAT1 test function
 * @param x:    Variables vector
 * @param f:    Objectives vector
 * @param g:    Constraints vector (no constraints, input: NULL)
 * @param xbin: Binary vector (no binary vector, input: NULL)
 */
void ZCAT11(double *x, double *f, double *g, int *xbin)
{
	int M = zcat_config.M;
	int n = zcat_config.N;
	int m = zcat_config.M - 1;
	double *y, *alpha, *beta;
	double *
	(*g_function)(double *y, int m, int n);
	void
	(*f_function)(double *F, double *y, int M);

	f_function = ZCAT_F[zcat11];
	g_function = (zcat_config.COMPLICATED_PS) ? ZCAT_G[zcat11] : g0;

	y = zcat_get_y(x, zcat_config.LB, zcat_config.UB, n); /* Normalization */
	alpha = zcat_get_alpha(y, M, f_function); /* Define Position */
	beta = zcat_get_beta(y, M, m, n, g_function); /* Define Distance*/
	zcat_mop_definition(f, alpha, beta, M); /* Assign fitness values */

	free(y);
	free(alpha);
	free(beta);
	return;
}

/**
 * ZCAT12 test function
 * @param x:    Variables vector
 * @param f:    Objectives vector
 * @param g:    Constraints vector (no constraints, input: NULL)
 * @param xbin: Binary vector (no binary vector, input: NULL)
 */
void ZCAT12(double *x, double *f, double *g, int *xbin)
{
	int M = zcat_config.M;
	int n = zcat_config.N;
	int m = zcat_config.M - 1;
	double *y, *alpha, *beta;
	double *
	(*g_function)(double *y, int m, int n);
	void
	(*f_function)(double *F, double *y, int M);

	f_function = ZCAT_F[zcat12];
	g_function = (zcat_config.COMPLICATED_PS) ? ZCAT_G[zcat12] : g0;

	y = zcat_get_y(x, zcat_config.LB, zcat_config.UB, n); /* Normalization */
	alpha = zcat_get_alpha(y, M, f_function); /* Define Position */
	beta = zcat_get_beta(y, M, m, n, g_function); /* Define Distance*/
	zcat_mop_definition(f, alpha, beta, M); /* Assign fitness values */

	free(y);
	free(alpha);
	free(beta);
	return;
}

/**
 * ZCAT13 test function
 * @param x:    Variables vector
 * @param f:    Objectives vector
 * @param g:    Constraints vector (no constraints, input: NULL)
 * @param xbin: Binary vector (no binary vector, input: NULL)
 */
void ZCAT13(double *x, double *f, double *g, int *xbin)
{
	int M = zcat_config.M;
	int n = zcat_config.N;
	int m = zcat_config.M - 1;
	double *y, *alpha, *beta;
	double *
	(*g_function)(double *y, int m, int n);
	void
	(*f_function)(double *F, double *y, int M);

	f_function = ZCAT_F[zcat13];
	g_function = (zcat_config.COMPLICATED_PS) ? ZCAT_G[zcat13] : g0;

	y = zcat_get_y(x, zcat_config.LB, zcat_config.UB, n); /* Normalization */
	alpha = zcat_get_alpha(y, M, f_function); /* Define Position */
	beta = zcat_get_beta(y, M, m, n, g_function); /* Define Distance*/
	zcat_mop_definition(f, alpha, beta, M); /* Assign fitness values */

	free(y);
	free(alpha);
	free(beta);
	return;
}

/**
 * ZCAT14 test function
 * @param x:    Variables vector
 * @param f:    Objectives vector
 * @param g:    Constraints vector (no constraints, input: NULL)
 * @param xbin: Binary vector (no binary vector, input: NULL)
 */
void ZCAT14(double *x, double *f, double *g, int *xbin)
{
	int M = zcat_config.M;
	int n = zcat_config.N;
	int m = 1;
	double *y, *alpha, *beta;
	double *
	(*g_function)(double *y, int m, int n);
	void
	(*f_function)(double *F, double *y, int M);

	f_function = ZCAT_F[zcat14];
	g_function = (zcat_config.COMPLICATED_PS) ? ZCAT_G[zcat14] : g0;

	y = zcat_get_y(x, zcat_config.LB, zcat_config.UB, n); /* Normalization */
	alpha = zcat_get_alpha(y, M, f_function); /* Define Position */
	beta = zcat_get_beta(y, M, m, n, g_function); /* Define Distance*/
	zcat_mop_definition(f, alpha, beta, M); /* Assign fitness values */

	free(y);
	free(alpha);
	free(beta);
	return;
}

/**
 * ZCAT15 test function
 * @param x:    Variables vector
 * @param f:    Objectives vector
 * @param g:    Constraints vector (no constraints, input: NULL)
 * @param xbin: Binary vector (no binary vector, input: NULL)
 */
void ZCAT15(double *x, double *f, double *g, int *xbin)
{
	int M = zcat_config.M;
	int n = zcat_config.N;
	int m = 1;
	double *y, *alpha, *beta;
	double *
	(*g_function)(double *y, int m, int n);
	void
	(*f_function)(double *F, double *y, int M);

	f_function = ZCAT_F[zcat15];
	g_function = (zcat_config.COMPLICATED_PS) ? ZCAT_G[zcat15] : g0;

	y = zcat_get_y(x, zcat_config.LB, zcat_config.UB, n); /* Normalization */
	alpha = zcat_get_alpha(y, M, f_function); /* Define Position */
	beta = zcat_get_beta(y, M, m, n, g_function); /* Define Distance*/
	zcat_mop_definition(f, alpha, beta, M); /* Assign fitness values */

	free(y);
	free(alpha);
	free(beta);
	return;
}

/**
 * ZCAT16 test function
 * @param x:    Variables vector
 * @param f:    Objectives vector
 * @param g:    Constraints vector (no constraints, input: NULL)
 * @param xbin: Binary vector (no binary vector, input: NULL)
 */
void ZCAT16(double *x, double *f, double *g, int *xbin)
{
	int M = zcat_config.M;
	int n = zcat_config.N;
	int m = 1;
	double *y, *alpha, *beta;
	double *
	(*g_function)(double *y, int m, int n);
	void
	(*f_function)(double *F, double *y, int M);

	f_function = ZCAT_F[zcat16];
	g_function = (zcat_config.COMPLICATED_PS) ? ZCAT_G[zcat16] : g0;

	y = zcat_get_y(x, zcat_config.LB, zcat_config.UB, n); /* Normalization */
	alpha = zcat_get_alpha(y, M, f_function); /* Define Position */
	beta = zcat_get_beta(y, M, m, n, g_function); /* Define Distance*/
	zcat_mop_definition(f, alpha, beta, M); /* Assign fitness values */

	free(y);
	free(alpha);
	free(beta);
	return;
}

/**
 * ZCAT17 test function
 * @param x:    Variables vector
 * @param f:    Objectives vector
 * @param g:    Constraints vector (no constraints, input: NULL)
 * @param xbin: Binary vector (no binary vector, input: NULL)
 */
void ZCAT17(double *x, double *f, double *g, int *xbin)
{
	int M = zcat_config.M;
	int n = zcat_config.N;
	int m = zcat_config.M - 1;
	double *y, *alpha, *beta;
	double *
	(*g_function)(double *y, int m, int n);
	void
	(*f_function)(double *F, double *y, int M);

	f_function = ZCAT_F[zcat17];
	g_function = (zcat_config.COMPLICATED_PS) ? ZCAT_G[zcat17] : g0;

	y = zcat_get_y(x, zcat_config.LB, zcat_config.UB, n); /* Normalization */
	alpha = zcat_get_alpha(y, M, f_function); /* Define Position */
	beta = zcat_get_beta(y, M, m, n, g_function); /* Define Distance*/
	zcat_mop_definition(f, alpha, beta, M); /* Assign fitness values */

	free(y);
	free(alpha);
	free(beta);
	return;
}

/**
 * ZCAT18 test function
 * @param x:    Variables vector
 * @param f:    Objectives vector
 * @param g:    Constraints vector (no constraints, input: NULL)
 * @param xbin: Binary vector (no binary vector, input: NULL)
 */
void ZCAT18(double *x, double *f, double *g, int *xbin)
{
	int M = zcat_config.M;
	int n = zcat_config.N;
	int m = zcat_config.M - 1;
	double *y, *alpha, *beta;
	double *
	(*g_function)(double *y, int m, int n);
	void
	(*f_function)(double *F, double *y, int M);

	f_function = ZCAT_F[zcat18];
	g_function = (zcat_config.COMPLICATED_PS) ? ZCAT_G[zcat18] : g0;

	y = zcat_get_y(x, zcat_config.LB, zcat_config.UB, n); /* Normalization */
	alpha = zcat_get_alpha(y, M, f_function); /* Define Position */
	beta = zcat_get_beta(y, M, m, n, g_function); /* Define Distance*/
	zcat_mop_definition(f, alpha, beta, M); /* Assign fitness values */

	free(y);
	free(alpha);
	free(beta);
	return;
}

/**
 * ZCAT19 test function
 * @param x:    Variables vector
 * @param f:    Objectives vector
 * @param g:    Constraints vector (no constraints, input: NULL)
 * @param xbin: Binary vector (no binary vector, input: NULL)
 */
void ZCAT19(double *x, double *f, double *g, int *xbin)
{
	int M = zcat_config.M;
	int n = zcat_config.N;
	int m;
	double *y, *alpha, *beta;
	double *
	(*g_function)(double *y, int m, int n);
	void
	(*f_function)(double *F, double *y, int M);

	f_function = ZCAT_F[zcat19];
	g_function = (zcat_config.COMPLICATED_PS) ? ZCAT_G[zcat19] : g0;

	y = zcat_get_y(x, zcat_config.LB, zcat_config.UB, n); /* Normalization */
	m = (zcat_value_in(y[0], 0.0, 0.2) || zcat_value_in(y[0], 0.4, 0.6)) ?
			1 : M - 1;
	alpha = zcat_get_alpha(y, M, f_function); /* Define Position */
	beta = zcat_get_beta(y, M, m, n, g_function); /* Define Distance*/
	zcat_mop_definition(f, alpha, beta, M); /* Assign fitness values */

	free(y);
	free(alpha);
	free(beta);
	return;
}

/**
 * ZCAT20 test function
 * @param x:    Variables vector
 * @param f:    Objectives vector
 * @param g:    Constraints vector (no constraints, input: NULL)
 * @param xbin: Binary vector (no binary vector, input: NULL)
 */
void ZCAT20(double *x, double *f, double *g, int *xbin)
{
	int M = zcat_config.M;
	int n = zcat_config.N;
	int m;
	double *y, *alpha, *beta;
	double *
	(*g_function)(double *y, int m, int n);
	void
	(*f_function)(double *F, double *y, int M);

	f_function = ZCAT_F[zcat20];
	g_function = (zcat_config.COMPLICATED_PS) ? ZCAT_G[zcat20] : g0;

	y = zcat_get_y(x, zcat_config.LB, zcat_config.UB, n); /* Normalization */
	m = (zcat_value_in(y[0], 0.1, 0.4) || zcat_value_in(y[0], 0.6, 0.9)) ?
			1 : zcat_config.M - 1;

	alpha = zcat_get_alpha(y, M, f_function); /* Define Position */
	beta = zcat_get_beta(y, M, m, n, g_function); /* Define Distance*/
	zcat_mop_definition(f, alpha, beta, M); /* Assign fitness values */

	free(y);
	free(alpha);
	free(beta);
}

/* ********************************************************************************************** */
/* ********************************************************************************************** */
/* ********************************************************************************************** */
/* ********************************************************************************************** */

void zcat_default_settings(int nobj)
{
	/** Benchmark Configuration **************************************************************** **/
	zcat_config.M = nobj;
	zcat_config.N = nobj * 10;
	zcat_config.LEVEL = 1;
	zcat_config.BIAS = ZCAT_FALSE;
	zcat_config.COMPLICATED_PS = ZCAT_TRUE;
	zcat_config.IMBALANCE = ZCAT_FALSE;

	zcat_config.LB = (double*) malloc(sizeof(double) * zcat_config.N);
	zcat_config.UB = (double*) malloc(sizeof(double) * zcat_config.N);
	zcat_set_bounds(zcat_config.LB, zcat_config.UB, zcat_config.N);
}

/**
 * Configuration function
 * @param nvar
 * @param nobj
 * @param Level
 * @param Bias
 * @param Imbalance
 * @param LB
 * @param UB
 */
void zcat_set(int nvar, int nobj, int Level, int Bias, int Complicated_PS,
		int Imbalance, double **LB, double **UB)
{
	int default_nvars = nobj * 10;

	assert(1 <= Level && Level <= 6);
	assert(Bias == ZCAT_TRUE || Bias == ZCAT_FALSE);
	assert(Complicated_PS == ZCAT_TRUE || Complicated_PS == ZCAT_FALSE);
	assert(nobj > 1);
	assert(nvar >= (nobj - 1) || nvar == -1); //TODO: It changed...

	/** Benchmark Configuration **************************************************************** **/
	zcat_config.LEVEL = Level;
	zcat_config.BIAS = Bias;
	zcat_config.COMPLICATED_PS = Complicated_PS;
	zcat_config.IMBALANCE = Imbalance;
	zcat_config.M = nobj;
	zcat_config.N = (nvar == -1) ? default_nvars : nvar;
	zcat_config.LB = (double*) malloc(sizeof(double) * zcat_config.N);
	zcat_config.UB = (double*) malloc(sizeof(double) * zcat_config.N);
	zcat_set_bounds(zcat_config.LB, zcat_config.UB, zcat_config.N);

	*LB = zcat_config.LB;
	*UB = zcat_config.UB;
	return;
}

/**
 * Unset structures used in ZCAT
 */
void zcat_unset()
{
	free(zcat_config.LB);
	free(zcat_config.UB);
	return;
}
