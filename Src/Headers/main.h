#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<omp.h>

#include"physics.h"

#include"mesh.h"
#include"metric.h"
#include"boundaries.h"

#include"limiters.h"

#include"const.h"
#include"macros.h"
#include"io.h"
#include"param.h"


// dummy
double x1, x2, x3;
int c1(int n, int i);
int c2(int n, int i, int j);
int c3(int n, int i, int j, int k);

/* Define pointers */
double *U, *U1, *U2, *U3;
double *Q, *Q1, *Q2, *Q3;

double start, delta;
double K;

/* Define number of grids */
int Nx1, Nx2, Nx3;

/* Define domain */
double x1max, x2max, x3max;
double x1min, x2min, x3min;

void Allocate_Array();

void New_Size();

void Initial();

int Boundaries(double *B);

void Integration();

int RK1D(double *u, double *q, double *q1, double *q2, int order);

int RK2D(double *u, double *q, double *q1, double *q2, int order);

int RK3D(double *u, double *q, double *q1, double *q2, int order);

int Flux1D(vec_ *v, lim_ *l, int *I);
                                   
int Flux2D(vec_ *v, lim_ *l, int *I);
                                   
int Flux3D(vec_ *v, lim_ *l, int *I);

int Hll(double *F, flx_ *f, int x);

int Hllc(double *F, flx_ *f, int x);

int AMATRIX1D(double *u, vec_ *v, int *I);
                                                            
int AMATRIX2D(double *u, vec_ *v, int *I);
                                              
int AMATRIX3D(double *u, vec_ *v, int *I);

int VECTOR(int pm, char flux, lim_ *l, flx_ *f, int *I);

int MxV(double *M, double *V, double *L);

void RoundGen(double *num);
