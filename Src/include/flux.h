/**
 * @file flux.h
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Numerical fluxes and solvers functions defined.
 */

typedef struct
{
	double uR[eq+1];
	double uL[eq+1];
	double qR[eq+1];
	double qL[eq+1];
	double fR[eq+1];
	double fL[eq+1];
	double lR;
	double lL;
}flx_;

void Numerical_Flux_F(double *F, int pm, int *I);

void Numerical_Flux_G(double *F, int pm, int *I);

void Numerical_Flux_H(double *F, int pm, int *I);

void Riemann_Solver(double *F, flx_ *f, int x);
