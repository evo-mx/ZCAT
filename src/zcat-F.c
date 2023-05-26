/*
 * zcat-F.c
 *
 *  Created on: Apr 1, 2022
 *      Author: Saul Zapotecas
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "zcat-tools.h"

/**
 * F1: Shape of the PF (as in the paper)
 * @param F:    The shape of the PF
 * @param y:    The first 'm' normalized decision variables
 * @param M:    The number of objectives
 */
void F1(double *F, double *y, int M)
{
	int i, j;

	F[0] = 1.0;
	for (i = 1; i <= M - 1; ++i)
	{
		F[0] *= sin(y[i - 1] * M_PI_2);
	}
	F[0] = zcat_fix_to_01(F[0]);
	assert(0 <= F[0] && F[0] <= 1.0);

	for (j = 2; j <= M - 1; ++j)
	{
		F[j - 1] = 1.0;
		for (i = 1; i <= M - j; ++i)
		{
			F[j - 1] *= sin(y[i - 1] * M_PI_2);
		}
		F[j - 1] *= cos(y[M - j + 1 - 1] * M_PI_2);
		F[j - 1] = zcat_fix_to_01(F[j - 1]);
		assert(0 <= F[j - 1] && F[j - 1] <= 1.0);
	}

	F[M - 1] = 1.0 - sin(y[0] * M_PI_2);
	F[M - 1] = zcat_fix_to_01(F[M - 1]);
	assert(0 <= F[M - 1] && F[M - 1] <= 1.0);
	return;
}

/**
 * F2: Shape of the PF (as in the paper)
 * @param F:    The shape of the PF
 * @param y:    The first 'm' normalized decision variables
 * @param M:    The number of objectives
 */
void F2(double *F, double *y, int M)
{
	int i, j;

	F[0] = 1.0;
	for (i = 1; i <= M - 1; ++i)
	{
		F[0] *= 1.0 - cos(y[i - 1] * M_PI_2);
	}
	assert(0 <= F[0] && F[0] <= 1.0);

	for (j = 2; j <= M - 1; ++j)
	{
		F[j - 1] = 1.0;
		for (i = 1; i <= M - j; ++i)
		{
			F[j - 1] *= 1.0 - cos(y[i - 1] * M_PI_2);
		}
		F[j - 1] *= 1.0 - sin(y[M - j + 1 - 1] * M_PI_2);
		assert(0 <= F[j - 1] && F[j - 1] <= 1.0);
	}

	F[M - 1] = 1.0 - sin(y[0] * M_PI_2);
	assert(0 <= F[M - 1] && F[M - 1] <= 1.0);
	return;
}

/**
 * F3: Shape of the PF (as in the paper)
 * @param F:    The shape of the PF
 * @param y:    The first 'm' normalized decision variables
 * @param M:    The number of objectives
 */
void F3(double *F, double *y, int M)
{
	int i, j;

	F[0] = 0.0;
	for (i = 1; i <= M - 1; ++i)
	{
		F[0] += y[i - 1];
	}
	F[0] = F[0] / (M - 1);
	assert(0 <= F[0] && F[0] <= 1.0);

	for (j = 2; j <= M - 1; ++j)
	{
		F[j - 1] = 0.0;
		for (i = 1; i <= M - j; ++i)
		{
			F[j - 1] += y[i - 1];
		}
		F[j - 1] += (1 - y[M - j + 1 - 1]);
		F[j - 1] = F[j - 1] / (M - j + 1);
		assert(0 <= F[j - 1] && F[j - 1] <= 1.0);
	}

	F[M - 1] = 1 - y[0];
	assert(0 <= F[M - 1] && F[M - 1] <= 1.0);
	return;
}

/**
 * F4: Shape of the PF (as in the paper)
 * @param F:    The shape of the PF
 * @param y:    The first 'm' normalized decision variables
 * @param M:    The number of objectives
 */
void F4(double *F, double *y, int M)
{
	int i, j;
	double sum;

	for (j = 1; j <= M - 1; ++j)
	{
		F[j - 1] = y[j - 1];
		assert(0 <= F[j - 1] && F[j - 1] <= 1.0);
	}

	sum = 0.0;
	for (i = 1; i <= M - 1; ++i)
	{
		sum += y[i - 1];
	}

	F[M - 1] = 1.0 - sum / (M - 1);
	assert(0 <= F[M - 1] && F[M - 1] <= 1.0);
	return;
}

/**
 * F5: Shape of the PF (as in the paper)
 * @param F:    The shape of the PF
 * @param y:    The first 'm' normalized decision variables
 * @param M:    The number of objectives
 */
void F5(double *F, double *y, int M)
{
	int i, j;
	double sum;

	for (j = 1; j <= M - 1; ++j)
	{
		F[j - 1] = y[j - 1];
		assert(0 <= F[j - 1] && F[j - 1] <= 1.0);
	}

	sum = 0.0;
	for (i = 1; i <= M - 1; ++i)
	{
		sum += 1.0 - y[i - 1];
	}

	F[M - 1] = (pow(exp(sum / (M - 1)), 8.0) - 1.0) / (pow(exp(1), 8.0) - 1.0);
	assert(0 <= F[M - 1] && F[M - 1] <= 1.0);
	return;
}

/**
 * F6: Shape of the PF (as in the paper)
 * @param F:    The shape of the PF
 * @param y:    The first 'm' normalized decision variables
 * @param M:    The number of objectives
 */
void F6(double *F, double *y, int M)
{
	int i, j;
	double mu;
	double k = 40.0;
	double r = 0.05;

	for (j = 1; j <= M - 1; ++j)
	{
		F[j - 1] = y[j - 1];
		assert(0 <= F[j - 1] && F[j - 1] <= 1.0);
	}

	mu = 0.0;
	for (i = 1; i <= M - 1; ++i)
	{
		mu += y[i - 1];
	}
	mu = mu / (M - 1);

	F[M - 1] = (pow(1 + exp(2 * k * mu - k), -1.0) - r * mu
			- pow(1 + exp(k), -1.0) + r)
			/ (pow(1 + exp(-k), -1.0) - pow(1 + exp(k), -1.0) + r);
	assert(0 <= F[M - 1] && F[M - 1] <= 1.0);
	return;
}

/**
 * F7: Shape of the PF (as in the paper)
 * @param F:    The shape of the PF
 * @param y:    The first 'm' normalized decision variables
 * @param M:    The number of objectives
 */
void F7(double *F, double *y, int M)
{
	int i, j;
	double sum;

	for (j = 1; j <= M - 1; ++j)
	{
		F[j - 1] = y[j - 1];
		assert(0 <= F[j - 1] && F[j - 1] <= 1.0);
	}

	sum = 0.0;
	for (i = 1; i <= M - 1; ++i)
	{
		sum += pow(0.5 - y[i - 1], 5);
	}
	sum = sum / (2 * (M - 1) * pow(0.5, 5));

	F[M - 1] = sum + 0.5;
	assert(0 <= F[M - 1] && F[M - 1] <= 1.0);
	return;
}

/**
 * F8: Shape of the PF (as in the paper)
 * @param F:    The shape of the PF
 * @param y:    The first 'm' normalized decision variables
 * @param M:    The number of objectives
 */
void F8(double *F, double *y, int M)
{
	int i, j;

	F[0] = 1.0;
	for (i = 1; i <= M - 1; ++i)
	{
		F[0] *= 1.0 - sin(y[i - 1] * M_PI_2);
	}
	F[0] = 1.0 - F[0];
	assert(0 <= F[0] && F[0] <= 1.0);

	for (j = 2; j <= M - 1; ++j)
	{
		F[j - 1] = 1.0;
		for (i = 1; i <= M - j; ++i)
		{
			F[j - 1] *= 1.0 - sin(y[i - 1] * M_PI_2);
		}
		F[j - 1] *= 1.0 - cos(y[M - j + 1 - 1] * M_PI_2);
		F[j - 1] = 1.0 - F[j - 1];
		assert(0 <= F[j - 1] && F[j - 1] <= 1.0);
	}

	F[M - 1] = cos(y[0] * M_PI_2);
	assert(0 <= F[M - 1] && F[M - 1] <= 1.0);
	return;
}

/**
 * F9: Shape of the PF (as in the paper)
 * @param F:    The shape of the PF
 * @param y:    The first 'm' normalized decision variables
 * @param M:    The number of objectives
 */
void F9(double *F, double *y, int M)
{
	int i, j;

	F[0] = 0.0;
	for (i = 1; i <= M - 1; ++i)
	{
		F[0] += sin(y[i - 1] * M_PI_2);
	}
	F[0] = F[0] / (M - 1);
	assert(0 <= F[0] && F[0] <= 1.0);

	for (j = 2; j <= M - 1; ++j)
	{
		F[j - 1] = 0.0;
		for (i = 1; i <= M - j; ++i)
		{
			F[j - 1] += sin(y[i - 1] * M_PI_2);
		}
		F[j - 1] += cos(y[M - j + 1 - 1] * M_PI_2);
		F[j - 1] = F[j - 1] / (M - j + 1);
		assert(0 <= F[j - 1] && F[j - 1] <= 1.0);
	}

	F[M - 1] = cos(y[0] * M_PI_2);
	assert(0 <= F[M - 1] && F[M - 1] <= 1.0);
	return;
}

/**
 * F10: Shape of the PF (as in the paper)
 * @param F:    The shape of the PF
 * @param y:    The first 'm' normalized decision variables
 * @param M:    The number of objectives
 */
void F10(double *F, double *y, int M)
{
	int j;
	double sum, r;

	sum = 0.0;
	r = 0.02;
	for (j = 1; j <= M - 1; ++j)
	{
		sum += 1 - y[j - 1];
		F[j - 1] = y[j - 1];
		assert(0 <= F[j - 1] && F[j - 1] <= 1.0);
	}

	F[M - 1] = (pow(r, -1) - pow(sum / (M - 1) + r, -1))
			/ (pow(r, -1) - pow(1 + r, -1));
	assert(0 <= F[M - 1] && F[M - 1] <= 1.0);
	return;
}

/**
 * F11: Shape of the PF (as in the paper)
 * @param F:    The shape of the PF
 * @param y:    The first 'm' normalized decision variables
 * @param M:    The number of objectives
 */
void F11(double *F, double *y, int M)
{
	int i, j;
	double k = 4.0;

	F[0] = 0.0;
	for (i = 1; i <= M - 1; ++i)
	{
		F[0] += y[i - 1];
	}
	F[0] = F[0] / (M - 1);
	assert(0 <= F[0] && F[0] <= 1.0);

	for (j = 2; j <= M - 1; ++j)
	{
		F[j - 1] = 0.0;
		for (i = 1; i <= M - j; ++i)
		{
			F[j - 1] += y[i - 1];
		}
		F[j - 1] += (1 - y[M - j + 1 - 1]);
		F[j - 1] = F[j - 1] / (M - j + 1);
		assert(0 <= F[j - 1] && F[j - 1] <= 1.0);
	}

	F[M - 1] = (cos((2 * k - 1) * y[0] * M_PI) + 2 * y[0] + 4 * k * (1 - y[0])
			- 1) / (4 * k);
	assert(0 <= F[M - 1] && F[M - 1] <= 1.0);
	return;
}

/**
 * F12: Shape of the PF (as in the paper)
 * @param F:    The shape of the PF
 * @param y:    The first 'm' normalized decision variables
 * @param M:    The number of objectives
 */
void F12(double *F, double *y, int M)
{
	int i, j;
	double k = 3.0;

	F[0] = 1.0;
	for (i = 1; i <= M - 1; ++i)
	{
		F[0] *= (1.0 - y[i - 1]);
	}
	F[0] = 1.0 - F[0];
	assert(0 <= F[0] && F[0] <= 1.0);

	for (j = 2; j <= M - 1; ++j)
	{
		F[j - 1] = 1.0;
		for (i = 1; i <= M - j; ++i)
		{
			F[j - 1] *= (1.0 - y[i - 1]);
		}
		F[j - 1] *= y[M - j + 1 - 1];
		F[j - 1] = 1.0 - F[j - 1];
		assert(0 <= F[j - 1] && F[j - 1] <= 1.0);
	}

	F[M - 1] = (cos((2 * k - 1) * y[0] * M_PI) + 2 * y[0] + 4 * k * (1 - y[0])
			- 1) / (4 * k);
	assert(0 <= F[M - 1] && F[M - 1] <= 1.0);
	return;
}

/**
 * F13: Shape of the PF (as in the paper)
 * @param F:    The shape of the PF
 * @param y:    The first 'm' normalized decision variables
 * @param M:    The number of objectives
 */
void F13(double *F, double *y, int M)
{
	int i, j;
	double k = 3.0;

	F[0] = 0.0;
	for (i = 1; i <= M - 1; ++i)
	{
		F[0] += sin(y[i - 1] * M_PI_2);
	}
	F[0] = 1.0 - F[0] / (M - 1.0);
	assert(0 <= F[0] && F[0] <= 1.0);

	for (j = 2; j <= M - 1; ++j)
	{
		F[j - 1] = 0.0;
		for (i = 1; i <= M - j; ++i)
		{
			F[j - 1] += sin(y[i - 1] * M_PI_2);
		}
		F[j - 1] += cos(y[M - j + 1 - 1] * M_PI_2);
		F[j - 1] = 1.0 - F[j - 1] / (M - j + 1);
		assert(0 <= F[j - 1] && F[j - 1] <= 1.0);
	}

	F[M - 1] = 1.0
			- (cos((2 * k - 1) * y[0] * M_PI) + 2 * y[0] + 4 * k * (1 - y[0])
					- 1) / (4.0 * k);
	assert(0 <= F[M - 1] && F[M - 1] <= 1.0);
	return;
}

/**
 * F14: Shape of the PF (as in the paper)
 * @param F:    The shape of the PF
 * @param y:    The first 'm' normalized decision variables
 * @param M:    The number of objectives
 */
void F14(double *F, double *y, int M)
{
	int j;

	F[0] = pow(sin(y[0] * M_PI_2), 2.0);
	assert(0 <= F[0] && F[0] <= 1.0);

	for (j = 2; j <= M - 2; ++j)
	{
		F[j - 1] = pow(sin(y[0] * M_PI_2), 2.0 + (j - 1.0) / (M - 2.0));
		assert(0 <= F[j - 1] && F[j - 1] <= 1.0);
	}

	if (M > 2)
	{
		F[M - 2] = 0.5 * (1 + sin(6 * y[0] * M_PI_2 - M_PI_2));
		assert(0 <= F[M - 2] && F[M - 2] <= 1.0);
	}

	F[M - 1] = cos(y[0] * M_PI_2);
	assert(0 <= F[M - 1] && F[M - 1] <= 1.0);
	return;
}

/**
 * F15: Shape of the PF (as in the paper)
 * @param F:    The shape of the PF
 * @param y:    The first 'm' normalized decision variables
 * @param M:    The number of objectives
 */
void F15(double *F, double *y, int M)
{
	int j;
	double k = 3.0;

	for (j = 1; j <= M - 1; ++j)
	{
		F[j - 1] = pow(y[0], 1.0 + (j - 1.0) / (4.0 * M));
		assert(0 <= F[j - 1] && F[j - 1] <= 1.0);
	}

	F[M - 1] = (cos((2 * k - 1) * y[0] * M_PI) + 2 * y[0] + 4 * k * (1 - y[0])
			- 1) / (4 * k);
	assert(0 <= F[M - 1] && F[M - 1] <= 1.0);
	return;
}

/**
 * F16: Shape of the PF (as in the paper)
 * @param F:    The shape of the PF
 * @param y:    The first 'm' normalized decision variables
 * @param M:    The number of objectives
 */
void F16(double *F, double *y, int M)
{
	int j;
	double k = 5;

	F[0] = sin(y[0] * M_PI_2);
	assert(0 <= F[0] && F[0] <= 1.0);

	for (j = 2; j <= M - 2; ++j)
	{
		F[j - 1] = pow(sin(y[0] * M_PI_2), 1.0 + (j - 1.0) / (M - 2.0));
		assert(0 <= F[j - 1] && F[j - 1] <= 1.0);
	}

	if (M > 2)
	{
		F[M - 2] = 0.5 * (1 + sin(10 * y[0] * M_PI_2 - M_PI_2));
		assert(0 <= F[M - 2] && F[M - 2] <= 1.0);
	}

	F[M - 1] = (cos((2 * k - 1) * y[0] * M_PI) + 2 * y[0] + 4 * k * (1 - y[0])
			- 1) / (4 * k);
	assert(0 <= F[M - 1] && F[M - 1] <= 1.0);
	return;
}

/**
 * F17: Shape of the PF (as in the paper)
 * @param F:    The shape of the PF
 * @param y:    The first 'm' normalized decision variables
 * @param M:    The number of objectives
 */
void F17(double *F, double *y, int M)
{
	int j;
	int wedge_flag;
	double sum;

	wedge_flag = zcat_forall_value_in(y, M - 1, 0.0, 0.5);

	sum = 0.0;
	for (j = 1; j <= M - 1; ++j)
	{
		if (wedge_flag)
		{
			F[j - 1] = y[0];
		}
		else
		{
			F[j - 1] = y[j - 1];
			sum += 1 - y[j - 1];
		}
		assert(0 <= F[j - 1] && F[j - 1] <= 1.0);
	}

	if (wedge_flag)
	{
		F[M - 1] = (pow(exp(1.0 - y[0]), 8.0) - 1.0)
				/ (pow(exp(1.0), 8.0) - 1.0);
	}
	else
	{
		F[M - 1] = (pow(exp(sum / (M - 1)), 8.0) - 1.0)
				/ (pow(exp(1.0), 8.0) - 1.0);
	}
	assert(0 <= F[M - 1] && F[M - 1] <= 1.0);
	return;
}

/**
 * F18: Shape of the PF (as in the paper)
 * @param F:    The shape of the PF
 * @param y:    The first 'm' normalized decision variables
 * @param M:    The number of objectives
 */
void F18(double *F, double *y, int M)
{
	int j;
	int wedge_flag, f1, f2;
	double sum;

	f1 = zcat_forall_value_in(y, M - 1, 0.0, 0.4);
	f2 = zcat_forall_value_in(y, M - 1, 0.6, 1.0);
	wedge_flag = (f1 == 1 || f2 == 1) ? 1 : 0;

	sum = 0.0;
	for (j = 1; j <= M - 1; ++j)
	{
		if (wedge_flag)
		{
			F[j - 1] = y[0];
		}
		else
		{
			F[j - 1] = y[j - 1];
			sum += pow(0.5 - y[j - 1], 5.0);
		}
		assert(0 <= F[j - 1] && F[j - 1] <= 1.0);
	}

	if (wedge_flag)
	{
		F[M - 1] = (pow(0.5 - y[0], 5.0) + pow(0.5, 5.0))
				/ (2.0 * pow(0.5, 5.0));
	}
	else
	{
		F[M - 1] = sum / (2.0 * (M - 1.0) * pow(0.5, 5.0)) + 0.5;
	}
	assert(0 <= F[M - 1] && F[M - 1] <= 1.0);
	return;
}

/**
 * F19: Shape of the PF (as in the paper) (mixed dimension 1 and M-1)
 * @param F:    The shape of the PF
 * @param y:    The first 'm' normalized decision variables
 * @param M:    The number of objectives
 */
void F19(double *F, double *y, int M)
{
	int j;
	double A = 5.0;
	double mu;

	int flag_deg;

	flag_deg =
			(zcat_value_in(y[0], 0.0, 0.2) || zcat_value_in(y[0], 0.4, 0.6)) ?
					1 : 0;
	mu = 0.0;
	for (j = 1; j <= M - 1; ++j)
	{
		mu += y[j - 1];
		F[j - 1] = (flag_deg == 1) ? y[0] : y[j - 1];
		F[j - 1] = zcat_fix_to_01(F[j - 1]);
		assert(0 <= F[j - 1] && F[j - 1] <= 1.0);
	}
	mu = (flag_deg == 1) ? y[0] : mu / (M - 1);

	F[M - 1] =
			(1.0 - mu - cos(2.0 * A * M_PI * mu + M_PI_2) / (2.0 * A * M_PI));
	F[M - 1] = zcat_fix_to_01(F[M - 1]);
	assert(0 <= F[M - 1] && F[M - 1] <= 1.0);
	return;
}

/**
 * F20: Shape of the PF (as in the paper)
 * @param F:    The shape of the PF
 * @param y:    The first 'm' normalized decision variables
 * @param M:    The number of objectives
 */
void F20(double *F, double *y, int M)
{
	int j;
	int deg_flag;
	double sum;

	deg_flag =
			(zcat_value_in(y[0], 0.1, 0.4) || zcat_value_in(y[0], 0.6, 0.9)) ?
					1 : 0;

	sum = 0.0;
	for (j = 1; j <= M - 1; ++j)
	{
		sum += pow(0.5 - y[j - 1], 5.0);
		F[j - 1] = (deg_flag == 1) ? y[0] : y[j - 1];
		assert(0 <= F[j - 1] && F[j - 1] <= 1.0);
	}

	if (deg_flag == 1)
	{
		F[M - 1] = (pow(0.5 - y[0], 5.0) + pow(0.5, 5.0))
				/ (2.0 * pow(0.5, 5.0));
	}
	else
	{
		F[M - 1] = sum / (2.0 * (M - 1.0) * pow(0.5, 5.0)) + 0.5;
	}
	assert(0 <= F[M - 1] && F[M - 1] <= 1.0);
	return;
}

