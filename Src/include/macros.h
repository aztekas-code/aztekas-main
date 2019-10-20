/**
 * @file macros.h
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Macros definitios for \a aztekas.
 *
 */

/* FUNCTIONS */
#define MIN(a,b) (((a)<(b))?(a):(b))                                            
#define MAX(a,b) (((a)>(b))?(a):(b)) 
#define Min(a,b,c,d)  (a < b ? a : b)  < (c < d ? c : d)  ? (a < b ? a : b)  : (c < d ? c : d)
#define Max(a,b,c,d)  (a > b ? a : b)  > (c > d ? c : d)  ? (a > b ? a : b)  : (c > d ? c : d)

/* PLUS AND MINUS FLAG */
#define MINUS 0
#define PLUS  1

/* LOGICAL */
#define FALSE 0
#define TRUE  1

/* PHYSICS */
#define HD    0
#define RHD   1

/* EOS */
#define IDEAL 0
#define DUST  1
#define STIFF 2

/* GRAVITY */
#define NONE   0
#define NEWTON 1

/* COORDINATES */
#define CARTESIAN   0
#define CYLINDRICAL 1
#define SPHERICAL   2

/* GRID */
#define UNIFORM   0
#define LOGMESH   1
#define TORTOISE  2   

/* FLUX RECONSTRUCTOR */
#define HLL  0
#define HLLC 1

/* INTEGRATION METHOD */
#define STANDARD 0
#define PVRS     1

/* RECONSTRUCTION */
#define GODUNOV  0
#define MINMOD   1
#define MC       2
#define SUPERBEE 3
#define WENO5    4

/* PRIMITIVE VARIABLES */
#define RHO       0
#define PRE       1
#define VX1       2
#define VX2       3
#define VX3       4

/* CONSERVATIVE VARIABLES */
#define DEN       0
#define ENE       1
#define MX1       2
#define MX2       3
#define MX3       4

/* SIMPLE ARRAY */
#define a(x,y)   a[(x)*(eq) + (y)]

#if DIM == 1
   #define  S1p(x)  grid.S1p[(x)]
   #define  S1m(x)  grid.S1m[(x)]

   #define  U(N,x)   U[(N)*(Nx1+1) + (x)]
   #define  U0(N,x)  U0[(N)*(Nx1+1) + (x)]
   #define  U1(N,x)  U1[(N)*(Nx1+1) + (x)]
   #define  U2(N,x)  U2[(N)*(Nx1+1) + (x)]

   #define  Q(N,x)   Q[(N)*(Nx1+1) + (x)]
   #define  Q0(N,x)  Q0[(N)*(Nx1+1) + (x)]
   #define  Q1(N,x)  Q1[(N)*(Nx1+1) + (x)]
   #define  Q2(N,x)  Q2[(N)*(Nx1+1) + (x)]

   #define  U1p(N,x)   U1p[(N)*(Nx1+1) + (x)]
   #define  U1m(N,x)   U1m[(N)*(Nx1+1) + (x)]

   #define  B(N,x)  B[(N)*(Nx1+1) + (x)]
   #define  u(N,x)  u[(N)*(Nx1+1) + (x)]
   #define  q(N,x)  q[(N)*(Nx1+1) + (x)]
   #define q1(N,x) q1[(N)*(Nx1+1) + (x)]
   #define q2(N,x) q2[(N)*(Nx1+1) + (x)]

#elif DIM == 2 || DIM == 4
   #define  S1p(x,y)  grid.S1p[(x)*(Nx2+1) + (y)]
   #define  S1m(x,y)  grid.S1m[(x)*(Nx2+1) + (y)]
   #define  S2p(x,y)  grid.S2p[(x)*(Nx2+1) + (y)]
   #define  S2m(x,y)  grid.S2m[(x)*(Nx2+1) + (y)]

   #define  U(N,x,y)  U[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define  U0(N,x,y) U0[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define  U1(N,x,y) U1[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define  U2(N,x,y) U2[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]

   #define  Q(N,x,y)  Q[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define  Q0(N,x,y) Q0[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define  Q1(N,x,y) Q1[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define  Q2(N,x,y) Q2[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]

   #define  U1p(N,x,y) U1p[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define  U1m(N,x,y) U1m[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define  U2p(N,x,y) U2p[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define  U2m(N,x,y) U2m[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]

   #define  B(N,x,y)  B[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define  u(N,x,y)  u[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define  q(N,x,y)  q[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define q1(N,x,y) q1[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define q2(N,x,y) q2[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]

#elif DIM == 3

   #define  U(N,x,y,z)  U[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]
   #define  Q(N,x,y,z)  Q[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]
   #define  B(N,x,y,z)  B[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]
   #define  u(N,x,y,z)  u[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]
   #define  q(N,x,y,z)  q[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]
   #define q1(N,x,y,z) q1[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]
   #define q2(N,x,y,z) q2[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]

#endif
