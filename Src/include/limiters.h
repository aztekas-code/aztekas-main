/**
 * @file limiters.h
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Reconstruction variables and functions definitions.
 *
 */

/*!
 * The structure \b lim_ contains vectors in which the reconstructed variables
 * of \f$U\f$ in each cell are stored.
 */
typedef struct
{
	double ux1p[2*eq];
	double ux1m[2*eq];
	double sx1[2*eq];
	double ux2p[2*eq];
	double ux2m[2*eq];
	double sx2[2*eq];
	double ux3p[2*eq];
	double ux3m[2*eq];
	double sx3[2*eq];
	double ux[2*eq];
   double U[eq+1];
}lim_;

void Primitive_Reconstruction();

double Limiter(double A, double B, int r);

double Godunov(double A, double B);

double Maxmod(double A, double B);

double Minmod(double A, double B);

double Mc(double A, double B);

double Superbee(double A, double B);

double Weno5(double v1, double v2, double v3, double v4, double v5);

int Reconst1D(double *u, lim_ *l, int *I);

int Reconst2D(double *u, lim_ *l, int *I);

int Reconst3D(double *u, lim_ *l, int *I);

