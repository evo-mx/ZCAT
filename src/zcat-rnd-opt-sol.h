/*
 * zcat-rnd-opt-sol.h
 *
 *  Created on: Apr 1, 2022
 *      Author: Saul Zapotecas
 */

#ifndef SRC_ZCAT_RND_OPT_SOL_H_
#define SRC_ZCAT_RND_OPT_SOL_H_

typedef struct
{
    double x1;
    double x2;
} SEGMENT;

double zcat_search_space(double y, int i);
void zcat_rnd_opt_sol(int mop, double *x, int nobj, int nreal);
SEGMENT *zcat_get_segments(int k);


#endif /* SRC_ZCAT_RND_OPT_SOL_H_ */
