/*
 * zca-rnd-opt-sol.c
 *
 *  Created on: Apr 1, 2022
 *      Author: Saul Zapotecas
 */

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "zcat-rnd-opt-sol.h"
#include "zcat-tools.h"
#include "zcat-benchmark.h"
#include "zcat-g.h"



/**
 * Mapping 'y' to the decision search space
 * @param y:    Normalize solution in [0,1]
 * @param i:    i-th decision variable
 * @return
 */
double zcat_search_space(double y, int i)
{
    double lb = -i * 0.5;
    double ub = i * 0.5;
    double x = y * (ub - lb) + lb;
    return x;
}

/**
 * Segments no dominated for disconnected MOPs
 * @param k: the number of disconnections (it depends on the MOP and the number of objectives)
 * @return
 */
SEGMENT *zcat_get_segments(int k)
{
    SEGMENT *segments;
    char fname[1024];
    FILE *fp = NULL;
    int i;

    segments = (SEGMENT*) malloc(sizeof(SEGMENT) * k);
    sprintf(fname, "seg/K%d.seg", k);

    fp = fopen(fname, "r");
    if (fp == NULL)
    {
        printf("ERROR: File %s no found\n", fname);
        exit(0);
    }

    for (i = 0; i < k; ++i)
    {
        if (fscanf(fp, "%lf", &segments[i].x1) == EOF)
        {
            printf("ERROR: Cannot read %s (get_segments).\n", fname);
            exit(0);
        }
        if (fscanf(fp, "%lf", &segments[i].x2) == EOF)
        {
            printf("ERROR: Cannot read %s (get_segments).\n", fname);
            exit(0);
        }
        //printf("[%e,%e]\n", segments[i].x1, segments[i].x2);
    }
    fclose(fp);
    return segments;
}

/**
 * Get a random optimal solution 'x' given a specific MOP
 * @param mop:  The reference of the mop
 * @param x:    The Pareto optimal solution
 * @param M:    The number of objectives
 * @param n:    The number of decision variables
 */
void zcat_rnd_opt_sol(int mop, double *x, int nobj, int nreal)
{
    int i, j;
    int m = -1;
    double y0, *g, *y;
    double *(*g_function)(double *, int, int);
    int k;
    SEGMENT * seg = NULL;

    g_function = ZCAT_G[mop];

    y0 = rnd_real(0.0, 1.0);
    m = nobj - 1;
    if ((zcat1 <= mop && mop <= zcat13) || zcat17 == mop || mop == zcat18)
    {
        m = nobj - 1;
    }
    else if (zcat14 <= mop && mop <= zcat16) // Degenerate PF
    {
        m = 1;
    }
    else if (mop == zcat19) // Hybrid
    {
        m = (zcat_value_in(y0, 0.0, 0.2) || zcat_value_in(y0, 0.4, 0.6)) ? 1 : nobj - 1;
    }
    else if (mop == zcat20) // Hybrid
    {
        m = (zcat_value_in(y0, 0.1, 0.4) || zcat_value_in(y0, 0.6, 0.9)) ? 1 : nobj - 1;
    }
    assert(m > 0);

    // For disconnected problems
    if ((zcat11 <= mop && mop <= zcat13) || mop == zcat15 || mop == zcat16)
    {
        /** Get the number of disconnections **/
        k = -1;
        if (mop == zcat11)
        {
            k = 4;
        }
        if (mop == zcat12 || mop == zcat13 || mop == zcat15)
        {
            k = 3;
        }
        if (mop == zcat16)
        {
            k = 5;
        }
        assert(k > 0);

        /** get  the range of each segment **/
        seg = zcat_get_segments(k);

        /** fix optimal solution **/
        k = rnd_int(0, k - 1); /** choose random segment **/
        y0 = rnd_real(seg[k].x1, seg[k].x2);
    }

    y = (double*) malloc(sizeof(double) * m);
    y[0] = y0;
    for (i = 2; i <= m; ++i)
    {
        y[i - 1] = rnd_real(0.0, 1.0);
    }

    g = g_function(y, m, nreal);

    /* Setting solution */
    for (i = 1; i <= m; ++i)
    {
        x[i - 1] = zcat_search_space(y[i - 1], i);
        assert((-i * 0.5) <= x[i - 1] && x[i - 1] <= (i * 0.5));
    }
    for (j = 0, i = m + 1; i <= nreal; ++i, ++j)
    {
        x[i - 1] = zcat_search_space(g[j], i);
        assert((-i * 0.5) <= x[i - 1] && x[i - 1] <= (i * 0.5));
    }

    free(seg);
    free(g);
    free(y);
    return;
}



