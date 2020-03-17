/**
 * @file main.h
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Main function, headers and variable declaration.
 */

#ifdef _OPENMP
   #include<omp.h>
   int MAX_NUM_THREADS;
   double start;
#else
   #include<time.h>
   clock_t start;
#endif

#include<unistd.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include"physics.h"

#include"boundaries.h"
#include"limiters.h"
#include"flux.h"

#include"const.h"
#include"macros.h"
#include"io.h"
#include"user_param.h"


// NAN
int CHECK_NAN;

/* Define pointers */
double *U, *U0, *U1, *U2, *U3;
double *Q, *Q0, *Q1, *Q2, *Q3;

double *U1p, *U1m;
double *U2p, *U2m;

double K;

/* Define number file */
int itprint;

/* Define freq. output dt and time */
double dtprint, tprint;


/* Define number of grids */
int Nx1, Nx2, Nx3;

/* Define domain */
double x1max, x2max, x3max;
double x1min, x2min, x3min;

void Allocate_Array();

void New_Size();

void Initial();

void Integration();

void Runge_Kutta(int order);

int MxV(double *M, double *V, double *L);

void RoundGen(double *num);
