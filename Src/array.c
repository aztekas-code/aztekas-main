/**
 * @file array.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Functions to simplify the index vector access.
 *
 * In this file we include three functions for passing the standard C
 * notation for a vector:
 * @code
 *    U[i*N_j*N_k + j*N_k + k]
 * @endcode
 * to a much simpler notation
 * @code
 *    U[(i,j,k)]
 * @endcode
 */

//Do not erase any of these libraries//
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include"main.h"
#include"param.h"

int c1(int n, int i)
{
   int c;

   c = n*(Nx1+1) + i;
   return c;
}

int c2(int n, int i, int j)
{
   int c;

   c = n*(Nx1+1)*(Nx2+1) + i*(Nx2+1) + j;
   return c;
}

int c3(int n, int i, int j, int k)
{
   int c;

   c = n*(Nx1+1)*(Nx2+1)*(Nx3+1) + i*(Nx2+1)*(Nx3+1) + j*(Nx3+1) + k;
   return c;
}
