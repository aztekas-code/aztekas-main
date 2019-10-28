/**
 * @file alloc.c
 * @author Alejandro Aguayo-Ortiz
 * @brief Essential allocation functions for \a aztekas.
 */

#include"main.h"

/*!
 * This function allocates the space in memory for all the vectors used in 
 * \a aztekas.
 */
void Allocate_Array()
{
   /**
    * Resize considering ghost cells
    */
   New_Size();

#if DIM == 1
   /**
    * Vectors for the X1 grid. X1p = X1(i+1/2) and X1m = X1(i-1/2)
    */
<<<<<<< HEAD
   grid.X1 = (double *)malloc((Nx1+1)*sizeof(double));
=======
   grid.X1  = (double *)malloc((Nx1+1)*sizeof(double));
>>>>>>> 12b3acd607466560c2bebb7b61677f23252c7907
   grid.X1p = (double *)malloc((Nx1+1)*sizeof(double));
   grid.X1m = (double *)malloc((Nx1+1)*sizeof(double));

   /**
    * Vectors for the interfaces surfaces S1(i+1/2)/dV and S1(i-1/2)/dV
    */
   grid.S1p = (double *)malloc((Nx1+1)*sizeof(double));
   grid.S1m = (double *)malloc((Nx1+1)*sizeof(double));

   /**
    * This is the solution vector in aztekas
    */
   U   = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double)); 
   Q   = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double)); 

   /**
    * Auxiliary vectors in time-stepping for the primitive variables U
    * and the conservative variables Q.
    */
   U0  = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double));
   U1  = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double));
   U2  = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double));
   Q0  = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double));
   Q1  = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double));
   Q2  = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double));

   /**
    * Vectors for the reconstructed primitive variables at the interfaces.
    */
   U1p  = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double));
   U1m  = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double));

#elif DIM == 2 || DIM == 4
   /**
    * Vectors for the X1 grid. X1p = X1(i+1/2) and X1m = X1(i-1/2)
    * Vectors for the X2 grid. X2p = X2(j+1/2) and X2m = X2(j-1/2)
    */
   grid.X1  = (double *)malloc((Nx1+1)*sizeof(double));
   grid.X1p = (double *)malloc((Nx1+1)*sizeof(double));
   grid.X1m = (double *)malloc((Nx1+1)*sizeof(double));
   grid.X2  = (double *)malloc((Nx2+1)*sizeof(double));
   grid.X2p = (double *)malloc((Nx2+1)*sizeof(double));
   grid.X2m = (double *)malloc((Nx2+1)*sizeof(double));

   /**
    * Vectors for the interfaces surfaces S1(i+1/2)/dV and S1(i-1/2)/dV
    * Vectors for the interfaces surfaces S2(j+1/2)/dV and S2(j-1/2)/dV
    */
   grid.S1p = (double *)malloc((Nx1+1)*(Nx2+1)*sizeof(double));
   grid.S1m = (double *)malloc((Nx1+1)*(Nx2+1)*sizeof(double));
   grid.S2p = (double *)malloc((Nx1+1)*(Nx2+1)*sizeof(double));
   grid.S2m = (double *)malloc((Nx1+1)*(Nx2+1)*sizeof(double));

   /**
    * This is the solution vector in aztekas
    */
   U   = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   Q   = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));

   /*
    * Auxiliary vectors in time-stepping for the primitive variables U
    * and the conservative variables Q.
    */
   U0  = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   U1  = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   U2  = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   Q0  = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   Q1  = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   Q2  = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));

   /**
    * Vectors for the reconstructed primitive variables at the interfaces.
    */
   U1p = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   U1m = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   U2p = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   U2m = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
#elif DIM == 3
   /**
    * Vectors for the X1 grid. X1p = X1(i+1/2) and X1m = X1(i-1/2)
    * Vectors for the X2 grid. X2p = X2(j+1/2) and X2m = X2(j-1/2)
    * Vectors for the X3 grid. X3p = X3(k+1/2) and X3m = X3(k-1/2)
    */
   grid.X1  = (double *)malloc((Nx1+1)*sizeof(double));
   grid.X1p = (double *)malloc((Nx1+1)*sizeof(double));
   grid.X1m = (double *)malloc((Nx1+1)*sizeof(double));
   grid.X2  = (double *)malloc((Nx2+1)*sizeof(double));
   grid.X2p = (double *)malloc((Nx2+1)*sizeof(double));
   grid.X2m = (double *)malloc((Nx2+1)*sizeof(double));
   grid.X3  = (double *)malloc((Nx3+1)*sizeof(double));
   grid.X3p = (double *)malloc((Nx3+1)*sizeof(double));
   grid.X3m = (double *)malloc((Nx3+1)*sizeof(double));

   /**
    * Vectors for the interfaces surfaces S1(i+1/2)/dV and S1(i-1/2)/dV
    * Vectors for the interfaces surfaces S2(j+1/2)/dV and S2(j-1/2)/dV
    * Vectors for the interfaces surfaces S3(k+1/2)/dV and S3(k-1/2)/dV
    */
   grid.S1p = (double *)malloc((Nx1+1)*(Nx2+1)*(Nx3+1)*sizeof(double));
   grid.S1m = (double *)malloc((Nx1+1)*(Nx2+1)*(Nx3+1)*sizeof(double));
   grid.S2p = (double *)malloc((Nx1+1)*(Nx2+1)*(Nx3+1)*sizeof(double));
   grid.S2m = (double *)malloc((Nx1+1)*(Nx2+1)*(Nx3+1)*sizeof(double));
   grid.S3p = (double *)malloc((Nx1+1)*(Nx2+1)*(Nx3+1)*sizeof(double));
   grid.S3m = (double *)malloc((Nx1+1)*(Nx2+1)*(Nx3+1)*sizeof(double));

   /**
    * This is the solution vector in aztekas
    */
   U   = (double *)malloc((Nx1+1)*(Nx2+1)*(Nx3+1)*(eq+1)*sizeof(double));

   /*
    * Auxiliary vectors in time-stepping for the primitive variables U
    * and the conservative variables Q.
    */
   U0  = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   U1  = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   U2  = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   Q0  = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   Q1  = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   Q2  = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));

   /**
    * Vectors for the reconstructed primitive variables at the interfaces.
    */
   U1p = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   U1m = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   U2p = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   U2m = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
#endif
}

/**
 * This function resize the length of all vectors considering the 
 * ghost cells. Now, if the mesh used to have 100 grid points, it will
 * now contain gc + 100 + gc grid points, where gc is the number of 
 * ghost cells (Default number = 3)
 */
void New_Size()
{
   Nx1 = Nx1 + 2*gc;
   Nx2 = Nx2 + 2*gc;
   Nx3 = Nx3 + 2*gc;
}
