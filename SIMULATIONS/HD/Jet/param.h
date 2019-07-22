// Define boundaries TRUE = 1 , FALSE = 0
// Set as cond_X whereas
// cond : {outflow,reflective,periodic,inflow}
// X : {x1max,x1min,x2max,x2min,x3max,x3min}

#define outflow_x1max      TRUE
#define reflective_x1min   TRUE
#define outflow_x2max      TRUE
#define reflective_x2min   TRUE

#define RECONST            MC
#define FLUX               HLL
#define GRID               UNIFORM

double n_jet, p_jet, vx1_jet, vx2_jet, vx3_jet;
double n_atm, p_atm, vx1_atm, vx2_atm, vx3_atm;
double r_jet, z_jet;
