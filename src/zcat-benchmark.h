/*
 * zcat-benchmark.h
 *
 *  Created on: Apr 1, 2022
 *      Author: Saul Zapotecas
 */

#ifndef BENCHMARCK_ZCAT_SZ_H_
#define BENCHMARCK_ZCAT_SZ_H_

enum
{
	zcat1,
	zcat2,
	zcat3,
	zcat4,
	zcat5,
	zcat6,
	zcat7,
	zcat8,
	zcat9,
	zcat10,
	zcat11,
	zcat12,
	zcat13,
	zcat14,
	zcat15,
	zcat16,
	zcat17,
	zcat18,
	zcat19,
	zcat20
};

extern char ZCAT_MOP_STR[20][7];
extern double *(*ZCAT_G[20])(double *y, int m, int n);
extern void (*ZCAT_F[20])(double *F, double *y, int M);
extern void (*ZCAT_MOP[20])(double *X, double *F, double *G, int *Xbin);

void ZCAT1(double *x, double *f, double *g, int *xbin);
void ZCAT2(double *x, double *f, double *g, int *xbin);
void ZCAT3(double *x, double *f, double *g, int *xbin);
void ZCAT4(double *x, double *f, double *g, int *xbin);
void ZCAT5(double *x, double *f, double *g, int *xbin);
void ZCAT6(double *x, double *f, double *g, int *xbin);
void ZCAT7(double *x, double *f, double *g, int *xbin);
void ZCAT8(double *x, double *f, double *g, int *xbin);
void ZCAT9(double *x, double *f, double *g, int *xbin);
void ZCAT10(double *x, double *f, double *g, int *xbin);
void ZCAT11(double *x, double *f, double *g, int *xbin);
void ZCAT12(double *x, double *f, double *g, int *xbin);
void ZCAT13(double *x, double *f, double *g, int *xbin);
void ZCAT14(double *x, double *f, double *g, int *xbin);
void ZCAT15(double *x, double *f, double *g, int *xbin);
void ZCAT16(double *x, double *f, double *g, int *xbin);
void ZCAT17(double *x, double *f, double *g, int *xbin);
void ZCAT18(double *x, double *f, double *g, int *xbin);
void ZCAT19(double *x, double *f, double *g, int *xbin);
void ZCAT20(double *x, double *f, double *g, int *xbin);

void zcat_set(int nvar, int nobj, int Level, int Bias, int Complicated_PS,
		int Imbalance, double **LB, double **UB);
void zcat_unset();

#endif /* BENCHMARCK_ZCAT_SZ_H_ */
