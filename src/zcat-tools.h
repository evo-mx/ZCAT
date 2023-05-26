/*
 * zcat-tools.h
 *
 *  Created on: Apr 1, 2022
 *      Author: Saul Zapotecas
 */

#ifndef SRC_ZCAT_TOOLS_H_
#define SRC_ZCAT_TOOLS_H_

double zcat_fix_to_01(double a);

int zcat_lq(double y, double z);
int zcat_eq(double y, double z);

int zcat_value_in(double y, double lb, double ub);
int zcat_forall_value_in(double *y, int m, double lb, double ub);

void print_double_vector(double *v, int size);
double rnd_real(double lb, double ub);
double rnd_perc();
int rnd_int(int lb, int ub);

#endif /* SRC_ZCAT_TOOLS_H_ */
