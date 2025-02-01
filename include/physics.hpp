/************************************************************************************
 *  @file     physics.hpp
 *
 *  @author   Alejandro Aguayo-Ortiz
 *
 *  @brief    Definition of functions, variables, and parameters used for
 *            all physical calculations. This C++ version takes advantage
 *            of language features such as standard containers for greater
 *            robustness and readability.
 *
 *  @details
 *  ---------------------------------------------------------------------------------
 *  This header declares the structures and functions related to the physics
 *  of the system. It includes the definition of two main structures:
 *    - rhs_: Contains arrays and values required to calculate different
 *            magnitudes (U, Q, F, etc.).
 *    - eos_: Defines essential thermodynamic values for the Equation of State.
 *
 *  In addition, it declares the functions that carry out conversions between
 *  primitive and conserved variables, flux calculations, source terms, metric
 *  operations, and calls to external (Fortran) routines for the Nuclear Equation of State.
 *
 *  IMPORTANT NOTE:
 *  ---------------------------------------------------------------------------------
 *  - All original comments have been kept, and additional comments have been
 *    added to clarify the functionality of certain parts of the code.
 *  - In C++, it is preferable to use `std::array` for fixed-size arrays rather
 *    than C-style arrays. This approach improves safety and readability.
 *  - Compatibility with `extern` external routines is maintained.
 *  - Ensure that the constant `eq` (or the macro `eq`) is a compile-time value
 *    so that `std::array<double, eq+1>` works properly.
 *
 ************************************************************************************/

#ifndef INCLUDE_PHYSICS_HPP_
#define INCLUDE_PHYSICS_HPP_

#include <array>

/**
 * @brief Structure that groups the necessary quantities for operations
 *        related to physical term and flux calculations.
 *
 * Each field is now a `std::array` of size eq+1. It also includes matrices and
 * auxiliary arrays for flux calculations in different directions (Fp, Fm, Gp, Gm, etc.),
 * and to store matrix A.
 *
 * It has been ported to C++ using `std::array` for better data management
 * safety and compile-time boundary checks.
 */
struct rhs_ {
    std::array<double, eq+1> U;   //!< Array of conserved variables
    std::array<double, eq+1> Q;   //!< Array of primitive variables
    std::array<double, eq+1> F;   //!< General flux array
    std::array<double, eq+1> L;   //!< Auxiliary array for term calculations
    std::array<double, (eq+1)*(eq+1)> A; //!< Matrix A of size (eq+1) x (eq+1)

    std::array<double, eq+1> Fp;  //!< Positive flux in the X direction
    std::array<double, eq+1> Fm;  //!< Negative flux in the X direction
    std::array<double, eq+1> Gp;  //!< Positive flux in the Y direction
    std::array<double, eq+1> Gm;  //!< Negative flux in the Y direction
    std::array<double, eq+1> Hp;  //!< Positive flux in the Z direction
    std::array<double, eq+1> Hm;  //!< Negative flux in the Z direction

    std::array<double, eq+1> S;   //!< Source terms
};

/**
 * @brief Structure for the Equation of State (EoS).
 *
 * It contains thermodynamic properties and their derivatives.
 * Each parameter is essential for describing the relationship between
 * density, pressure, internal energy, entropy, and other thermodynamic
 * quantities.
 */
struct eos_ {
    double rho;     //!< Density
    double p;       //!< Pressure
    double e;       //!< Internal energy
    double s;       //!< Entropy
    double temp;    //!< Temperature
    double cs;      //!< Speed of sound
    double dpt;     //!< Partial derivative of pressure w.r.t. temperature
    double det;     //!< Partial derivative of energy w.r.t. temperature
    double h;       //!< Enthalpy
    double dhdrho;  //!< Partial derivative of enthalpy w.r.t. density
    double dhdp;    //!< Partial derivative of enthalpy w.r.t. pressure
};

/**
 * @brief Adiabatic index (Gamma), used in polytropic equations of state.
 *        Must be declared extern so it can be used in other files.
 */
extern double K;


/****************************************************************************************
 *                            Function Declarations
 ***************************************************************************************/

/**
 * @brief Auxiliary function to compute matrix A or part of it.
 *
 * @param a   Pointer to an array of coefficients (output).
 * @param uu  Pointer to an array of conserved or primitive variables (input).
 *
 * @return An integer indicating success (0) or failure.
 */
int funct_A(double *a, double *uu);

/**
 * @brief Converts conserved variables to primitive variables.
 *
 * @param q Pointer to array of primitive variables (output).
 * @param u Pointer to array of conserved variables (input).
 *
 * @return An integer indicating success (0) or failure.
 */
int Cons2Prim(double *q, double *u);

/**
 * @brief Converts all primitive variables to conserved variables for the entire domain.
 *
 * @param u Pointer to array of conserved variables (output).
 * @param q Pointer to array of primitive variables (input).
 */
void Prim2Cons_All(double *u, double *q);

/**
 * @brief Converts primitive variables to conserved variables in a specific local region.
 *
 * @param q           Pointer to array of primitive variables (input).
 * @param u           Pointer to array of conserved variables (output).
 * @param local_grid  Pointer to the local grid or gauge for the conversion.
 */
void Prim2Cons(double *q, double *u, gauge_ *local_grid);

/**
 * @brief Computes flux F in the X direction from primitive variables.
 *
 * @param f           Pointer to flux array (output).
 * @param v           Pointer to additional primitive variables (input).
 * @param u           Pointer to conserved variables (input).
 * @param local_grid  Pointer to the local grid or gauge.
 */
void Prim2FluxF(double *f, double *v, double *u, gauge_ *local_grid);

/**
 * @brief Computes flux G in the Y direction from primitive variables.
 *
 * @param f           Pointer to flux array (output).
 * @param v           Pointer to additional primitive variables (input).
 * @param u           Pointer to conserved variables (input).
 * @param local_grid  Pointer to the local grid or gauge.
 */
void Prim2FluxG(double *f, double *v, double *u, gauge_ *local_grid);

/**
 * @brief Computes flux H in the Z direction from primitive variables.
 *
 * @param f           Pointer to flux array (output).
 * @param v           Pointer to additional primitive variables (input).
 * @param u           Pointer to conserved variables (input).
 * @param local_grid  Pointer to the local grid or gauge.
 */
void Prim2FluxH(double *f, double *v, double *u, gauge_ *local_grid);

/**
 * @brief Converts primitive variables to source terms.
 *
 * @param s Pointer to an array of source terms (output).
 * @param I Pointer to an integer or index array (input).
 */
void Prim2Sources(double *s, int *I);

/**
 * @brief Calculates source terms from conserved variables.
 *
 * @param s           Pointer to an array of source terms (output).
 * @param u           Pointer to an array of conserved variables (input).
 * @param local_grid  Pointer to the local grid or gauge.
 */
void Source_Terms(double *s, double *u, gauge_ *local_grid);

/**
 * @brief Calculates user-defined source terms to customize the physics.
 *
 * @param s           Pointer to an array of source terms (output).
 * @param u           Pointer to an array of conserved variables (input).
 * @param local_grid  Pointer to the local grid or gauge.
 */
void User_Source_Terms(double *s, double *u, gauge_ *local_grid);

/**
 * @brief Builds or updates the matrix A that plays a role in the system of equations.
 *
 * @param a           Pointer to matrix A (output).
 * @param u           Pointer to conserved variables (input).
 * @param local_grid  Local grid or gauge (passed by value or reference).
 */
void Matrix_A(double *a, double *u, gauge_ local_grid);

/**
 * @brief Calculates the equation of state from conserved variables.
 *
 * @param e           Pointer to the eos_ structure (output).
 * @param u           Pointer to conserved variables (input).
 * @param local_grid  Pointer to the local grid or gauge.
 */
void EoS(eos_ *e, double *u, gauge_ *local_grid);

/**
 * @brief Calculates derivatives with respect to temperature for the nuclear EoS.
 *
 * @param e  Pointer to the eos_ structure (output).
 * @param u  Pointer to conserved variables (input).
 */
void EoS_DT(eos_ *e, double *u);

/**
 * @brief Obtains the metric components used in the local gauge.
 *
 * @param local_grid Pointer to the gauge_ structure where the components are stored.
 */
void Get_Metric_Components(gauge_ *local_grid);

/**
 * @brief Calculates gauge derivatives, possibly for curvature or shift terms.
 *
 * @param der         Pointer to the der_gauge_ structure where derivatives are stored.
 * @param local_grid  Pointer to the local grid or gauge.
 */
void Gauge_Derivatives(der_gauge_ *der, gauge_ *local_grid);

/**
 * @brief Computes the scalar contraction of rank 1 (for example, a dot product
 *        in a curved space).
 *
 * @param scalar      Pointer to the location where the resulting scalar is stored.
 * @param cov         Pointer to a covariant array.
 * @param con         Pointer to a contravariant array.
 */
void Scalar_Contraction_Range1(double *scalar, double *cov, double *con);

/**
 * @brief Raises the index of a rank-1 vector using the inverse metric.
 *
 * @param con         Pointer to the contravariant vector (output).
 * @param cov         Pointer to the covariant vector (input).
 * @param local_grid  Metric information in the local grid.
 */
void Raise_Index_Range1(double *con, double *cov, gauge_ *local_grid);

/**
 * @brief Lowers the index of a rank-1 vector using the metric.
 *
 * @param cov         Pointer to the covariant vector (output).
 * @param con         Pointer to the contravariant vector (input).
 * @param local_grid  Metric information in the local grid.
 */
void Low_Index_Range1(double *cov, double *con, gauge_ *local_grid);

/**
 * @brief Lowers indices of a rank-2 object (for example, a matrix or tensor).
 *
 * @param diag        Covariant matrix (output).
 * @param con         Contravariant matrix (input).
 * @param local_grid  Metric information in the local grid.
 */
void Low_Index_Range2(double **diag, double **con, gauge_ *local_grid);


/****************************************************************************************
 *                            External (Fortran) Routines
 ***************************************************************************************/
/**
 * @brief Fortran routine for the nuclear EoS calculation with respect to T (temperature).
 *
 * @param var       Input variables array.
 * @param xxMass    Mass or mass fraction array.
 * @param AA        Mass number.
 * @param ZZ        Atomic number.
 * @param term_var  Output array for EoS-related terms.
 */
extern void nad_eos_dt_(double var[], double xxMass[], double AA[], double ZZ[], double term_var[]);

/**
 * @brief Fortran routine for the nuclear EoS calculation with respect to P (pressure).
 *
 * @param var       Input variables array.
 * @param xxMass    Mass or mass fraction array.
 * @param AA        Mass number.
 * @param ZZ        Atomic number.
 * @param term_var  Output array for EoS-related terms.
 */
extern void nad_eos_dp_(double var[], double xxMass[], double AA[], double ZZ[], double term_var[]);

/**
 * @brief Fortran routine for the nuclear EoS calculation with respect to E (energy).
 *
 * @param var       Input variables array.
 * @param xxMass    Mass or mass fraction array.
 * @param AA        Mass number.
 * @param ZZ        Atomic number.
 * @param term_var  Output array for EoS-related terms.
 */
extern void nad_eos_de_(double var[], double xxMass[], double AA[], double ZZ[], double term_var[]);

/**
 * @brief Fortran routine for the calculation of nuclear observables (or other magnitudes).
 *
 * @param term_var  Output array for EoS-related terms.
 */
extern void nados_(double term_var[]);

#endif // INCLUDE_PHYSICS_HPP_
