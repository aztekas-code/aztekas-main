/************************************************************************************
 *  @file     flux.hpp
 *
 * .@author   Alejandro Aguayo-Ortiz
 *  @brief    Numerical fluxes and solver function definitions for aztekas.
 *
 *  @details
 *  ---------------------------------------------------------------------------------
 *  This header provides:
 *    - A structure (`flx_`) for storing left and right states (uR, uL, etc.) used in
 *      Riemann solvers and flux computations.
 *    - Function prototypes for various numerical flux functions and the Riemann solver.
 *    - Uses `std::array<double, kEq + 1>` for fixed-size arrays, where kEq maps to `eq`.
 *
 *  IMPORTANT NOTES:
 *  ---------------------------------------------------------------------------------
 *  - We assume `eq` is a compile-time constant (e.g., `#define eq 5` or `constexpr int eq = 5;`).
 *    We define `kEq` for clarity and consistency.
 *  - If `eq` is not known at compile time, you must switch to a dynamic container like
 *    `std::vector<double>` and refactor accordingly.
 *
 ************************************************************************************/

#ifndef INCLUDE_FLUX_HPP_
#define INCLUDE_FLUX_HPP_

#include <array>

// Map your compile-time constant eq to a named constant kEq
// e.g. #define eq 5  =>  enum { kEq = eq };
enum { kEq = eq };

/**
 * @brief Structure used for flux calculations, storing left/right states and
 *        associated fluxes.
 *
 * - `uR`, `uL`: Right and left conserved variable states.
 * - `qR`, `qL`: Right and left primitive variable states.
 * - `fR`, `fL`: Right and left fluxes.
 * - `lR`, `lL`: Possibly characteristic speeds or wave speeds from the Riemann solver.
 */
struct flx_ {
    std::array<double, kEq + 1> uR;  //!< Right state (conserved vars)
    std::array<double, kEq + 1> uL;  //!< Left state (conserved vars)

    std::array<double, kEq + 1> qR;  //!< Right state (primitive vars)
    std::array<double, kEq + 1> qL;  //!< Left state (primitive vars)

    std::array<double, kEq + 1> fR;  //!< Right flux
    std::array<double, kEq + 1> fL;  //!< Left flux

    double lR; //!< Possibly the maximum wave speed to the right
    double lL; //!< Possibly the maximum wave speed to the left
};


/****************************************************************************************
 *                           Function Declarations
 ***************************************************************************************/

/**
 * @brief Computes the numerical flux in the X1 direction (often referred to as F-flux).
 *
 * @param F   Pointer to the output flux array.
 * @param pm  Sign or direction indicator (plus/minus).
 * @param I   Pointer to an integer array of indices (i, j, k).
 */
void Numerical_Flux_F(double* F, int pm, int* I);

/**
 * @brief Computes the numerical flux in the X2 direction (often referred to as G-flux).
 *
 * @param F   Pointer to the output flux array.
 * @param pm  Sign or direction indicator (plus/minus).
 * @param I   Pointer to an integer array of indices (i, j, k).
 */
void Numerical_Flux_G(double* F, int pm, int* I);

/**
 * @brief Computes the numerical flux in the X3 direction (often referred to as H-flux).
 *
 * @param F   Pointer to the output flux array.
 * @param pm  Sign or direction indicator (plus/minus).
 * @param I   Pointer to an integer array of indices (i, j, k).
 */
void Numerical_Flux_H(double* F, int pm, int* I);

/**
 * @brief Executes the Riemann solver given the flux structure (`flx_`) and an index.
 *
 * @param F   Pointer to the output flux array.
 * @param f   Pointer to a `flx_` struct containing left/right states.
 * @param x   Integer specifying the coordinate index or direction.
 */
void Riemann_Solver(double* F, flx_* f, int x);

#endif // INCLUDE_FLUX_HPP_
