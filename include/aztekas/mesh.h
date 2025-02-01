#ifndef INCLUDE_AZTEKAS_MESH_H_
#define INCLUDE_AZTEKAS_MESH_H_

#include <vector>

namespace Mesh {

// Define number of grids
extern int Nx1, Nx2, Nx3;

// Define domain limits
extern double x1max, x2max, x3max;
extern double x1min, x2min, x3min;

// Grid class
class Grid {
public:
  double time;

  std::vector<double> X1, X1p, X1m;
  std::vector<double> X2, X2p, X2m;
  std::vector<double> X3, X3p, X3m;

  std::vector<double> S1p, S1m;
  std::vector<double> S2p, S2m;
  std::vector<double> S3p, S3m;

  Grid(int size1D, int size2D, int size3D) {
    X1.resize(size1D);
    X1p.resize(size1D);
    X1m.resize(size1D);
    X2.resize(size2D);
    X2p.resize(size2D);
    X2m.resize(size2D);
    X3.resize(size3D);
    X3p.resize(size3D);
    X3m.resize(size3D);

    S1p.resize(size1D);
    S1m.resize(size1D);
    S2p.resize(size2D);
    S2m.resize(size2D);
    S3p.resize(size3D);
    S3m.resize(size3D);
  }
};

extern Grid grid;

// Gauge class
class Gauge {
public:
  double x[4];
  double lapse;
  double beta_con[3];
  double gamma_con[3][3];
  double dety;
};

// Derivative of Gauge class
class DerivativeGauge {
public:
  double dlapse[3];
  double dbeta[3][3];
  double dgam[3][3][3];
};

// Define additional parameters
extern double dx1, dx2, dx3;
extern double dt;
extern double tmax, cou;

// Function prototypes
int Mesh();
void Surface_Volume();
double TimeStep();

} // namespace Mesh

#endif // INCLUDE_AZTEKAS_MESH_H_
