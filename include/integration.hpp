/************************************************************************************
 *  @file     integration.hpp
 *
 *  @author   Alejandro Aguayo-Ortiz
 *  @brief    Integration functions, structures, and global variables for aztekas.
 *
 *  @details
 *  ---------------------------------------------------------------------------------
 *  - Defines the Runge-Kutta (rk_) structure with fields for intermediate steps
 *    in integration (u0, k1, etc.).
 *  - Declares global variables for solution arrays (U, Q, etc.) using std::vector.
 *  - Declares function prototypes for solving equation systems, hyperbolic integration,
 *    ODE integration, and more.
 *
 *  IMPORTANT NOTES:
 *  ---------------------------------------------------------------------------------
 *  - The actual memory allocation/initialization for the global std::vectors (U, U0, ...)
 *    should be done in a corresponding .cpp file (e.g., integration.cpp).
 *  - If `rk_order` or `itprint` are used across multiple translation units, they should
 *    also be defined in exactly one .cpp file without `extern`.
 *
 ************************************************************************************/

#ifndef INCLUDE_INTEGRATION_HPP_
#define INCLUDE_INTEGRATION_HPP_

#include <vector>  // For std::vector

/**
 * @brief A structure for Runge-Kutta integration. It contains intermediate
 *        integration values (u0, k1, k2, etc.), a step size 'h', and a function
 *        value 'f' if needed.
 */
struct rk_ {
    double u0;
    double u1;
    double u2;
    double u3;
    double u4;

    double k1;
    double k2;
    double k3;
    double k4;

    double h;  //!< Step size for integration
    double f;  //!< Function value or auxiliary variable
};

/**
 * @brief Order of the Runge-Kutta method (e.g., 2 for RK2, 4 for RK4).
 *        Declared extern; define in a .cpp file.
 */
extern int rk_order;

/*-----------------------------------------------------------------------------------
 * Global solution arrays
 * Using std::vector for safer memory management.
 * You can resize these vectors at runtime according to your needs (e.g., Nx1*Nx2*Nx3).
 *----------------------------------------------------------------------------------*/
extern std::vector<double> U, U0, U1, U2, U3;
extern std::vector<double> Q, Q0, Q1, Q2, Q3;

extern std::vector<double> U1p, U1m;
extern std::vector<double> U2p, U2m;

/**
 * @brief The current iteration for printing (or checkpointing).
 */
extern int itprint;

/**
 * @brief Frequency of output in terms of time (dtprint) and the current time for output (tprint).
 */
extern double dtprint, tprint;


/*-----------------------------------------------------------------------------------
 * Function Prototypes
 *----------------------------------------------------------------------------------*/

/**
 * @brief Solves the system of equations (possibly spatial derivatives, source terms, etc.).
 */
void Equation_System_Solver();

/**
 * @brief Performs hyperbolic integration (e.g., explicit time stepping of PDEs).
 */
void Hyperbolic_Integration();

/**
 * @brief Generic ODE integration routine that might handle advanced/time-dependent ODEs.
 */
void ODE_Integration();

/**
 * @brief A Runge-Kutta integrator that takes a pointer to rk_ structure and the chosen order.
 *
 * @param rk    Pointer to a structure holding intermediate Runge-Kutta values.
 * @param order The order of Runge-Kutta (rk_order).
 */
void Runge_Kutta(rk_* rk, int order);

/**
 * @brief Method of Lines approach to PDE integration, using the specified order (order).
 *
 * @param order The order of accuracy for the integration method (likely 2, 3, or 4).
 */
void Method_of_Lines(int order);

/**
 * @brief Matrix-vector product routine, possibly for linear algebra steps (e.g., M*x = L).
 *
 * @param M Pointer to a matrix in row-major form.
 * @param V Pointer to a vector.
 * @param L Pointer to the result vector (output).
 *
 * @return An integer indicating success/failure.
 */
int MxV(double* M, double* V, double* L);

#endif  // INCLUDE_INTEGRATION_HPP_
