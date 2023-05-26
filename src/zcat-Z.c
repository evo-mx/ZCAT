/*
 * zcat-Z.c
 *
 *  Created on: Apr 1, 2022
 *      Author: Saul Zapotecas
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

/**
 * Z1 function
 * @param J:        Subset of decision variables to consider from YII - g(YI) (cf the paper)
 * @param Jsize:    The size of J
 * @return
 */
double Z1(double *J, int Jsize)
{
    double Z;
    int i;

    assert(Jsize > 0);
    Z = 0.0;
    for (i = 0; i < Jsize; ++i)
    {
        Z += J[i] * J[i];
    }
    Z = (10.0 / Jsize) * Z;
    assert(0 <= Z && Z <= 10.0);
    return (Z);
}

/**
 * Z2 function
 * @param J:        Subset of decision variables to consider from YII - g(YI) (cf the paper)
 * @param Jsize:    The size of J
 * @return
 */
double Z2(double *J, int Jsize)
{
    double Z;
    int i;

    assert(Jsize > 0);
    Z = -DBL_MAX;
    for (i = 0; i < Jsize; ++i)
    {
        Z = fmax(Z, fabs(J[i]));
    }
    Z = 10.0 * Z;
    assert(0 <= Z && Z <= 10.0);
    return (Z);
}

/**
 * Z3 function
 * @param J:        Subset of decision variables to consider from YII - g(YI) (cf the paper)
 * @param Jsize:    The size of J
 * @return
 */
double Z3(double *J, int Jsize)
{
    int i;
    double Z, k;

    assert(Jsize > 0);
    k = 5.0;
    Z = 0.0;
    for (i = 0; i < Jsize; ++i)
    {
        Z += (pow(J[i], 2.0) - cos((2.0 * k - 1) * M_PI * J[i]) + 1.0) / 3.0;
    }
    Z = (10.0 / Jsize) * Z;
    assert(0 <= Z && Z <= 10.0);
    return (Z);
}

/**
 * Z4 function
 * @param J:        Subset of decision variables to consider from YII - g(YI) (cf the paper)
 * @param Jsize:    The size of J
 * @return
 */
double Z4(double *J, int Jsize)
{
    int i;
    double Z, k, pow1, pow2;

    assert(Jsize > 0);
    k = 5.0;
    Z = 0.0;
    pow1 = -DBL_MAX;
    pow2 = 0.0;
    for (i = 0; i < Jsize; ++i)
    {
        pow1 = fmax(pow1, fabs(J[i]));
        pow2 += 0.5 * (cos((2.0 * k - 1) * M_PI * J[i]) + 1.0);
    }
    Z = (10.0 / (2.0 * exp(1) - 2.0)) * (exp(pow(pow1, 0.5)) - exp(pow2 / Jsize) - 1.0 + exp(1));
    assert(0 <= Z && Z <= 10.0);
    return (Z);
}

/**
 * Z5 function
 * @param J:        Subset of decision variables to consider from YII - g(YI) (cf the paper)
 * @param Jsize:    The size of J
 * @return
 */
double Z5(double *J, int Jsize)
{
    int i;
    double Z;

    assert(Jsize > 0);
    Z = 0.0;
    for (i = 0; i < Jsize; ++i)
    {
        Z += pow(fabs(J[i]), 0.002);
    }
    Z = -0.7 * Z3(J, Jsize) + (10.0 / Jsize) * Z;
    assert(0 <= Z && Z <= 10.0);
    return (Z);
}

/**
 * Z6 function
 * @param J:        Subset of decision variables to consider from YII - g(YI) (cf the paper)
 * @param Jsize:    The size of J
 * @return
 */
double Z6(double *J, int Jsize)
{
    int i;
    double Z;

    assert(Jsize > 0);
    Z = 0.0;
    for (i = 0; i < Jsize; ++i)
    {
        Z += fabs(J[i]);
    }
    Z = -0.7 * Z4(J, Jsize) + 10.0 * pow(Z / Jsize, 0.002);
    assert(0 <= Z && Z <= 10.0);
    return (Z);
}

/**
 * Bias transformation
 * @param z:    value to bias
 * @return
 */
double Zbias(double z)
{
    double w;
    double gamma = 0.05;

    w = pow(fabs(z), gamma);
    assert(0.0 <= w && w <= 1.0);
    return w;
}
