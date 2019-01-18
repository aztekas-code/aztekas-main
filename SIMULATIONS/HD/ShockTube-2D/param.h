// Define boundaries TRUE = 1 , FALSE = 0
// Set as cond_X whereas
// cond : {outflow,reflective,periodic,inflow}
// X : {x1max,x1min,x2max,x2min,x3max,x3min} for outflow, reflective and inflow
// X : {x1,x2,x3} for periodic
#define outflow_x1max 1
#define outflow_x1min 1
#define outflow_x2max 1
#define outflow_x2min 1

// Define interface position
// interface = 0 --> Along X
// interface = 1 --> Along Y
// interface = 2 --> Along a diagonal x + y - 1 = 0
#define interface 0

#define limiter 'W'
#define riemann 2
