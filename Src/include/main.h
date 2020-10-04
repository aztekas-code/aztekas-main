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

#include"mesh.h"
#include"physics.h"

#include"initial.h"
#include"integration.h"
#include"boundaries.h"
#include"limiters.h"
#include"flux.h"

#include"const.h"
#include"macros.h"
#include"io.h"
#include"user_param.h"
