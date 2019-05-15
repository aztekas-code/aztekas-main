void RIEMANN();
void KH();
void JET();
void SPH_ACC();
void WIND();
void INIT_CUSTOM();

void OUTFLOW(double *B);
void PERIODIC(double *B);
void REFLECTIVE(double *B);
void JET_LAUNCH(double *B);
void IN_OUT_BOUND(double *B);
void WIND_BOUND(double *B);
void BOUND_CUSTOM();

void RESTART();
void RESTART_BIN();
void EXTFORCE(double *a, double *uu);
