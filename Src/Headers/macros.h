/* FUNCTIONS */
#define MIN(a,b) (((a)<(b))?(a):(b))                                            
#define MAX(a,b) (((a)>(b))?(a):(b)) 

/* LOGICAL */
#define TRUE  1
#define FALSE 0

/* PHYSICS */
#define HD    0
#define RHD   1

/* EOS */
#define IDEAL 0

/* COORDINATES */
#define CARTESIAN   0
#define CYLINDRICAL 1
#define SPHERICAL   2

/* GRID */
#define UNIFORM 0
#define LOGMESH 1

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

/* SIMPLE ARRAY */
#if DIM == 1
   #define  S1p(x)  grid.S1p[(x)]
   #define  S1m(x)  grid.S1m[(x)]

   #define  U(N,x)  U[(N)*(Nx1+1) + (x)]
   #define  Q(N,x)  Q[(N)*(Nx1+1) + (x)]
   #define  B(N,x)  B[(N)*(Nx1+1) + (x)]
   #define  u(N,x)  u[(N)*(Nx1+1) + (x)]
   #define  a(N,x)  a[(N)*(Nx1+1) + (x)]
   #define uu(N,x) uu[(N)*(Nx1+1) + (x)]
   #define  q(N,x)  q[(N)*(Nx1+1) + (x)]
   #define q1(N,x) q1[(N)*(Nx1+1) + (x)]
   #define q2(N,x) q2[(N)*(Nx1+1) + (x)]

#elif DIM == 2 || DIM == 4
   #define  S1p(x,y)  grid.S1p[(x)*(Nx2+1) + (y)]
   #define  S1m(x,y)  grid.S1m[(x)*(Nx2+1) + (y)]
   #define  S2p(x,y)  grid.S2p[(x)*(Nx2+1) + (y)]
   #define  S2m(x,y)  grid.S2m[(x)*(Nx2+1) + (y)]

   #define  U(N,x,y)  U[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define  Q(N,x,y)  Q[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define  B(N,x,y)  B[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define  u(N,x,y)  u[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define  a(N,x,y)  a[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define uu(N,x,y) uu[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define  q(N,x,y)  q[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define q1(N,x,y) q1[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]
   #define q2(N,x,y) q2[(N)*(Nx1+1)*(Nx2+1) + (x)*(Nx2+1) + (y)]

#elif DIM == 3

   #define  U(N,x,y,z)  U[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]
   #define  Q(N,x,y,z)  Q[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]
   #define  B(N,x,y,z)  B[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]
   #define  u(N,x,y,z)  u[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]
   #define  a(N,x,y,z)  a[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]
   #define uu(N,x,y,z) uu[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]
   #define  q(N,x,y,z)  q[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]
   #define q1(N,x,y,z) q1[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]
   #define q2(N,x,y,z) q2[(N)*(Nx1+1)*(Nx2+1)*(Nx3+1) + (x)*(Nx2+1)*(Nx3+1) + (y)*(Nx3+1) + (z)]

#endif

