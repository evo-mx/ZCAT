/*
 * zcat-tools.c
 *
 *  Created on: Apr 1, 2022
 *      Author: Saul Zapotecas
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

/**
 * Fix the value 'a' to be between 0 and 1 (numerical precision issues)
 * @param a:    The value to fix
 * @return:     The fixed value
 */
double zcat_fix_to_01(double a)
{
	double epsilon = DBL_EPSILON;
	double min = 0.0;
	double max = 1.0;

	double min_epsilon = min - epsilon;
	double max_epsilon = max + epsilon;

	if (a <= min && a >= min_epsilon)
	{
		return min;
	}
	else if (a >= max && a <= max_epsilon)
	{
		return max;
	}
	else
	{
		return a;
	}
}

/**
 * For precision issues.
 * This function verify whether y<=z
 */
int zcat_lq(double y, double z)
{
	if (y < z)
	{
		return 1;
	}
	if (fabs(z - y) < DBL_EPSILON)
	{
		return 1;
	}
	return 0;
}

/**
 * For precision issues.
 * This function verify whether y==z
 */
int zcat_eq(double y, double z)
{
	if (fabs(z - y) < DBL_EPSILON)
	{
		return 1;
	}
	return 0;
}


/**
 * Verify if value 'y' is into 'lb' and 'ub': lb<= y <= ub
 * @param y:    A given value
 * @param lb:   Low bound
 * @param ub:   Up bound
 * @return: 1: lb<= y <= ub;
 *          0: otherwise
 */
int zcat_value_in(double y, double lb, double ub)
{
	return (zcat_lq(lb, y) && zcat_lq(y, ub)) ? 1 : 0;
}

/**
 * Verify if all the values 'yi' are into  'lb' and 'ub': lb<= yi <= ub (i=0,..,m-1)
 * @param y:    A given vector
 * @param m:    The size of y
 * @param lb:   Low bound
 * @param ub:   Up bound
 * @return: 1: lb<= yi <= ub (for all i=1,..,m-1);
 *          0: otherwise
 */
int zcat_forall_value_in(double *y, int m, double lb, double ub)
{
	int i;
	for (i = 0; i < m; ++i)
	{
		if (zcat_value_in(y[i], lb, ub) == 0)
		{
			return 0;
		}
	}
	return 1;
}

/**
 * Display double vector
 * @param v
 * @param size
 */
void print_double_vector(double *v, int size)
{
	int i;
	for (i = 0; i < size; ++i)
	{
		printf("%lf ", v[i]);
	}
	printf("\n");
	return;
}

/**
 * Real number between lb and ub (it includes the bounds)
 * @param lb: the low bound
 * @param ub: the up bound
 * @return
 */
double rnd_real(double lb, double ub)
{
	// rand() returns a pseudo-random integer between 0 and RAND_MAX
	double rnd = ((double) rand() / (double) (RAND_MAX));
	assert(lb < ub);
	rnd = (rnd * (ub - lb) + lb);
	assert(lb <= rnd && rnd <= ub);
	return rnd;
}

double rnd_perc()
{
	return rnd_real(0.0, 1.0);
}

/**
 * Integer number between lb and ub (it includes the bounds)
 * @param lb: the low bound
 * @param ub: the up bound
 * @return
 */
int rnd_int(int lb, int ub)
{
	// rand() returns a pseudo-random integer between 0 and RAND_MAX
	double rnd = ((double) rand() / (double) (RAND_MAX));
	int r = (int) (lb + (rnd * (ub - lb + 1)));
	assert(lb <= r && r <= ub);
	return r;
}

