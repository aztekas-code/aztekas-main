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
}lim_;

double Limiter(double A, double B, int r);

double Godunov(double A, double B);

double Maxmod(double A, double B);

double Minmod(double A, double B);

double Mc(double A, double B);

double Superbee(double A, double B);

double Weno5(double v1, double v2, double v3, double v4, double v5);

int Reconst1D(double *u, int r, lim_ *l, int *I);

int Reconst2D(double *u, int r, lim_ *l, int *I);

int Reconst3D(double *u, int r, lim_ *l, int *I);

