#ifndef INCLUDE_AZTEKAS_PHYSICS_H_
#define INCLUDE_AZTEKAS_PHYSICS_H_

#include "aztekas/mesh.h"
#include <array>

namespace Physics {

constexpr int EqSize = eq + 1;

class RHS {
public:
  std::array<double, EqSize> U, Q, F, L, Fp, Fm, Gp, Gm, Hp, Hm, S;
  std::array<double, EqSize * EqSize> A;
};

class EOS {
public:
  double rho, p, e, s, temp, cs, dpt, det, h, dhdrho, dhdp;
};

// Adiabatic index
extern double K;

// Function prototypes
int funct_A(double *a, double *uu);
int Cons2Prim(double *q, double *u);

void Prim2Cons_All(double *u, double *q);
void Prim2Cons(double *q, double *u, Mesh::Gauge &local_grid);
void Prim2FluxF(double *f, double *v, double *u, Mesh::Gauge &local_grid);
void Prim2FluxG(double *f, double *v, double *u, Mesh::Gauge &local_grid);
void Prim2FluxH(double *f, double *v, double *u, Mesh::Gauge &local_grid);

void Prim2Sources(double *s, int *I);
void Source_Terms(double *s, double *u, Mesh::Gauge &local_grid);
void User_Source_Terms(double *s, double *u, Mesh::Gauge &local_grid);

void Matrix_A(double *a, double *u, Mesh::Gauge &local_grid);
void EoS(EOS &e, double *u, Mesh::Gauge &local_grid);
void EoS_DT(EOS &e, double *u);

void Get_Metric_Components(Mesh::Gauge &local_grid);
void Gauge_Derivatives(Mesh::DerivativeGauge &der, Mesh::Gauge &local_grid);
void Scalar_Contraction_Range1(double *scalar, double *cov, double *con);
void Raise_Index_Range1(double *con, double *cov, Mesh::Gauge &local_grid);
void Low_Index_Range1(double *cov, double *con, Mesh::Gauge &local_grid);
void Low_Index_Range2(double **diag, double **con, Mesh::Gauge &local_grid);

// External functions
extern void nad_eos_dt_(double var[], double xxMass[], double AA[], double ZZ[],
                        double term_var[]);
extern void nad_eos_dp_(double var[], double xxMass[], double AA[], double ZZ[],
                        double term_var[]);
extern void nad_eos_de_(double var[], double xxMass[], double AA[], double ZZ[],
                        double term_var[]);
extern void nados_(double term_var[]);

} // namespace Physics

#endif // INCLUDE_AZTEKAS_PHYSICS_H_
