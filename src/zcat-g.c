/*
 * zcat-g.c
 *
 *  Created on: Apr 1, 2022
 *      Author: Saul Zapotecas
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/**
 * Angle for each dimension of PS
 * @param j:    Coordinate of PS
 * @param m:    Dimension of PS
 * @param n:    Number of decision variables of the MOP
 * @return
 */
double Thetaj(int j, int m, int n)
{
    double Tj;
    assert(1 <= j && j <= n - m);
    Tj = 2.0 * M_PI * (j - 1.0) / (n - m);
    return Tj;
}

/**
 * g0: Shape of the PS (exactly as DTLZx, WFGx, i.e., it is a piecewise linear PS)
 * @param y:    The first 'm' normalized decision variables
 * @param m:    The dimension of PS
 * @param n:    The number of decision variables
 * @return
 */
double *g0(double *y, int m, int n)
{
    int j;
    double *g;

    g = (double*) malloc(sizeof(double) * (n - m));
    for (j = 1; j <= n - m; ++j)
    {
        g[j - 1] = 0.2210;
        assert(0 <= g[j - 1] && g[j - 1] <= 1.0);
    }
    return (g);
}

/**
 * g1: Shape of the PS (as in the paper)
 * @param y:    The first 'm' normalized decision variables
 * @param m:    The dimension of PS
 * @param n:    The number of decision variables
 * @return
 */
double *g1(double *y, int m, int n)
{
    int i, j;
    double sum, *g;

    g = (double*) malloc(sizeof(double) * (n - m));
    for (j = 1; j <= n - m; ++j)
    {
        sum = 0.0;
        for (i = 1; i <= m; ++i)
        {
            sum += sin(1.5 * M_PI * y[i - 1] + Thetaj(j, m, n));
        }
        g[j - 1] = sum / (2.0 * m) + 0.5;
        assert(0 <= g[j - 1] && g[j - 1] <= 1.0);
    }
    return (g);
}

/**
 * g2: Shape of the PS (as in the paper)
 * @param y:    The first 'm' normalized decision variables
 * @param m:    The dimension of PS
 * @param n:    The number of decision variables
 * @return
 */
double *g2(double *y, int m, int n)
{
    int i, j;
    double sum, *g;

    g = (double*) malloc(sizeof(double) * (n - m));
    for (j = 1; j <= n - m; ++j)
    {
        sum = 0.0;
        for (i = 1; i <= m; ++i)
        {
            sum += pow(y[i - 1], 2.0) * sin(4.5 * M_PI * y[i - 1] + Thetaj(j, m, n));
        }
        g[j - 1] = sum / (2.0 * m) + 0.5;
        assert(0 <= g[j - 1] && g[j - 1] <= 1.0);
    }
    return (g);
}

/**
 * g3: Shape of the PS (as in the paper)
 * @param y:    The first 'm' normalized decision variables
 * @param m:    The dimension of PS
 * @param n:    The number of decision variables
 * @return
 */
double *g3(double *y, int m, int n)
{
    int i, j;
    double sum, *g;

    g = (double*) malloc(sizeof(double) * (n - m));
    for (j = 1; j <= n - m; ++j)
    {
        sum = 0.0;
        for (i = 1; i <= m; ++i)
        {
            sum += pow(cos(M_PI * y[i - 1] + Thetaj(j, m, n)), 2.0);
        }
        g[j - 1] = sum / m;
        assert(0 <= g[j - 1] && g[j - 1] <= 1.0);
    }
    return (g);
}

/**
 * g4: Shape of the PS (as in the paper)
 * @param y:    The first 'm' normalized decision variables
 * @param m:    The dimension of PS
 * @param n:    The number of decision variables
 * @return
 */
double *g4(double *y, int m, int n)
{
    int i, j;
    double mu, *g;

    g = (double*) malloc(sizeof(double) * (n - m));
    for (j = 1; j <= n - m; ++j)
    {
        mu = 0.0;
        for (i = 1; i <= m; ++i)
        {
            mu += y[i - 1];
        }
        mu = mu / m;
        g[j - 1] = (mu / 2.0) * cos(4.0 * M_PI * mu + Thetaj(j, m, n)) + 0.5;
        assert(0 <= g[j - 1] && g[j - 1] <= 1.0);
    }
    return (g);
}

/**
 * g5: Shape of the PS (as in the paper)
 * @param y:    The first 'm' normalized decision variables
 * @param m:    The dimension of PS
 * @param n:    The number of decision variables
 * @return
 */
double *g5(double *y, int m, int n)
{
    int i, j;
    double sum, *g;

    g = (double*) malloc(sizeof(double) * (n - m));
    for (j = 1; j <= n - m; ++j)
    {
        sum = 0.0;
        for (i = 1; i <= m; ++i)
        {
            sum += pow(sin(2.0 * M_PI * y[i - 1] - 1 + Thetaj(j, m, n)), 3.0);
        }
        g[j - 1] = sum / (2.0 * m) + 0.5;
        assert(0 <= g[j - 1] && g[j - 1] <= 1.0);
    }
    return (g);
}

/**
 * g6: Shape of the PS (as in the paper)
 * @param y:    The first 'm' normalized decision variables
 * @param m:    The dimension of PS
 * @param n:    The number of decision variables
 * @return
 */
double *g6(double *y, int m, int n)
{
    int i, j;
    double s1, s2, *g;

    g = (double*) malloc(sizeof(double) * (n - m));
    for (j = 1; j <= n - m; ++j)
    {
        s1 = 0.0;
        s2 = 0.0;
        for (i = 1; i <= m; ++i)
        {
            s1 += pow(y[i - 1], 2.0);
            s2 += pow(cos(11.0 * M_PI * y[i - 1] + Thetaj(j, m, n)), 3.0);
        }
        s1 /= m;
        s2 /= m;
        g[j - 1] = (-10.0 * exp((-2.0 / 5.0) * sqrt(s1)) - exp(s2) + 10.0 + exp(1.0))
                / (-10.0 * exp(-2.0 / 5.0) - pow(exp(1.0), -1) + 10.0 + exp(1.0));
        assert(0 <= g[j - 1] && g[j - 1] <= 1.0);
    }
    return (g);
}

/**
 * g7: Shape of the PS (as in the paper)
 * @param y:    The first 'm' normalized decision variables
 * @param m:    The dimension of PS
 * @param n:    The number of decision variables
 * @return
 */
double *g7(double *y, int m, int n)
{
    int i, j;
    double mu, *g;

    g = (double*) malloc(sizeof(double) * (n - m));
    for (j = 1; j <= n - m; ++j)
    {
        mu = 0.0;
        for (i = 1; i <= m; ++i)
        {
            mu += y[i - 1];
        }
        mu /= m;
        g[j - 1] = (mu + exp(sin(7.0 * M_PI * mu - M_PI_2 + Thetaj(j, m, n))) - exp(-1.0))
                / (1.0 + exp(1) - exp(-1));
        assert(0 <= g[j - 1] && g[j - 1] <= 1.0);
    }
    return (g);
}

/**
 * g8: Shape of the PS (as in the paper)
 * @param y:    The first 'm' normalized decision variables
 * @param m:    The dimension of PS
 * @param n:    The number of decision variables
 * @return
 */
double *g8(double *y, int m, int n)
{
    int i, j;
    double sum, *g;

    g = (double*) malloc(sizeof(double) * (n - m));
    for (j = 1; j <= n - m; ++j)
    {
        sum = 0.0;
        for (i = 1; i <= m; ++i)
        {
            sum += fabs(sin(2.5 * M_PI * (y[i - 1] - 0.5) + Thetaj(j, m, n)));
        }
        g[j - 1] = sum / m;
        assert(0 <= g[j - 1] && g[j - 1] <= 1.0);
    }
    return (g);
}

/**
 * g9: Shape of the PS (as in the paper)
 * @param y:    The first 'm' normalized decision variables
 * @param m:    The dimension of PS
 * @param n:    The number of decision variables
 * @return
 */
double *g9(double *y, int m, int n)
{
    int i, j;
    double mu, sum, *g;

    g = (double*) malloc(sizeof(double) * (n - m));
    for (j = 1; j <= n - m; ++j)
    {
        sum = 0.0;
        mu = 0.0;
        for (i = 1; i <= m; ++i)
        {
            sum += fabs(sin(2.5 * M_PI * y[i - 1] - M_PI_2 + Thetaj(j, m, n)));
            mu += y[i - 1];
        }
        mu /= m;
        g[j - 1] = mu / 2.0 - sum / (2.0 * m) + 0.5;
        assert(0 <= g[j - 1] && g[j - 1] <= 1.0);
    }
    return (g);
}

/**
 * g10: Shape of the PS (as in the paper)
 * @param y:    The first 'm' normalized decision variables
 * @param m:    The dimension of PS
 * @param n:    The number of decision variables
 * @return
 */
double *g10(double *y, int m, int n)
{
    int i, j;
    double sum, *g;

    g = (double*) malloc(sizeof(double) * (n - m));
    for (j = 1; j <= n - m; ++j)
    {
        sum = 0.0;
        for (i = 1; i <= m; ++i)
        {
            sum += sin((4.0 * y[i - 1] - 2.0) * M_PI + Thetaj(j, m, n));
        }
        g[j - 1] = pow(sum, 3.0) / (2.0 * pow(m, 3.0)) + 0.5;
        assert(0 <= g[j - 1] && g[j - 1] <= 1.0);
    }
    return (g);
}

