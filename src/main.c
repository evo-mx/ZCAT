/*
 * main.c
 *
 *  Created on: Apr 1, 2022
 *      Author: Saul Zapotecas
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <time.h>

#include "zcat-tools.h"
#include "zcat-benchmark.h"
#include "zcat-rnd-opt-sol.h"


void print_file_header(FILE *file, char *mop, int nobj, int nvar)
{
	int j;
	fprintf(file, "# %s: %d objectives, %d decision variables\n# ", mop, nobj,
			nvar);
	for (j = 0; j < nobj; ++j)
		fprintf(file, "<obj%d> ", j + 1);
	for (j = 0; j < nvar; ++j)
		fprintf(file, "<var%d> ", j + 1);
	fprintf(file, "\n");
	return;
}

void print_file_solution(FILE *file, double *x, double *f, int nobj, int nvar)
{
	int j;
	for (j = 0; j < nobj; ++j)
		fprintf(file, "%e ", f[j]);
	for (j = 0; j < nvar; ++j)
		fprintf(file, "%e ", x[j]);
	fprintf(file, "\n");
	return;
}

/* main function */
int main()
{
	int i;
	int nobj = 3; /* number of objectives */
	int nvar = 10 * nobj; /* standard decision variables */

	int Level = 1; /* Level of the problem {1,..,6} */
	int Bias_flag = 0; /* Bias flag {0,1}  (True: 1, False: 0)*/
	int Complicated_PS_flag = 1; /* Complicated PS flag {0,1} (True: 1, False: 0) */
	int Imbalance_flag = 0; /* Imbalance flag {0,1}  (True: 1, False: 0)*/

	double *LB = NULL, *UB = NULL; /* Pointer to the bounds of the problem */
	double *x, *f;	/* Decision variables (x) and objectives (f)*/

	FILE *pf;
	int max_solutions;
	int mop;
	char fname[1024];

	void (*test_problem)(double*, double*, double*, int *); /* test problem */

	/* seed for random numbers */
	srand(time(NULL));

	/* Configuring ZCAT benchmark. LB and UB are pointers to the bounds in the ZCAT structures */
	zcat_set(nvar, nobj, Level, Bias_flag, Complicated_PS_flag, Imbalance_flag, &LB, &UB);

	assert(LB!=NULL && UB!=NULL);
	x = (double*) malloc(sizeof(double) * nvar);
	f = (double*) malloc(sizeof(double) * nobj);

	printf("\nLB: ");
	print_double_vector(LB, nvar);
	printf("UB: ");
	print_double_vector(UB, nvar);

	/* Example 1: Generating and evaluating a single random solution */
	for (i = 0; i < nvar; ++i) /* Generating random solution */
	{
		x[i] = rnd_real(LB[i], UB[i]);
	}
	ZCAT1(x, f, NULL, NULL); /* Evaluating random solution */
	printf("===================\n");
	printf("Example 0\n");
	printf("===================\n");
	printf("Decision variables: \n");
	print_double_vector(x, nvar);
	printf("Objective values: \n");
	print_double_vector(f, nobj);
	printf("===================\n");

	max_solutions = 600; /* Maximum number of optimal solutions */

	/* Example 2: How to generate random optimal solutions and how to use the ZCATx functions to
	 * evaluate any solution. The following lines generate 'max_solutions' optimal solutions and
	 * write them into the file: Example-ZCAT1-<nobj>objs.opt (see the resulting file) */
	sprintf(fname, "optimal-solutions/Example0-%s-%dobjs.opt",
			ZCAT_MOP_STR[zcat1], nobj);
	pf = fopen(fname, "w");
	print_file_header(pf, ZCAT_MOP_STR[zcat1], nobj, nvar);
	for (i = 0; i < max_solutions; ++i)
	{
		zcat_rnd_opt_sol(zcat1, x, nobj, nvar); /* Generating optimal solution */
		ZCAT1(x, f, NULL, NULL); /* Here, ZCAT1 evaluates the generated optimal solution.
		 Since ZCAT does not consider constraints and binary decision variables,
		 arguments 3 and 4 have to be NULL.
		This is the way in which any solutions should be evaluated */

		print_file_solution(pf, x, f, nobj, nvar); /* print solution in the file */
	}
	fclose(pf);

	/* Example 3: Automatic generation of 'max_solutions' optimal solutions for all the problems.
	 * The following lines generate 'max_solutions' optimal solutions for each problem and write
	 * them into the file: ZCAT<1,..,20>-<nobj>objs.opt (see the resulting files)
	 * */
	for (mop = zcat1; mop <= zcat20; mop++)
	{
		sprintf(fname, "optimal-solutions/%s-%dobjs.opt", ZCAT_MOP_STR[mop],
				nobj);
		pf = fopen(fname, "w");
		print_file_header(pf, ZCAT_MOP_STR[mop], nobj, nvar);
		for (i = 0; i < max_solutions; ++i)
		{
			zcat_rnd_opt_sol(mop, x, nobj, nvar); /* Generating optimal solution */
			test_problem = ZCAT_MOP[mop];
			test_problem(x, f, NULL, NULL);
			print_file_solution(pf, x, f, nobj, nvar); /* print solution in the file */
		}
		fclose(pf);
	}
	zcat_unset(); /* Necessary to free the bounds in ZCAT structures.
	 It could be replaced with: free(zcat_config.LB) and free(zcat_config.UB) */
	free(x);
	free(f);
	return 0;
}

