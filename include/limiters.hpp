/************************************************************************************
 *  @file     limiters.hpp
 *
 * .@author   Alejandro Aguayo-Ortiz
 *  @brief    Definitions of reconstruction variables and function prototypes
 *            related to limiters, using modern C++ containers.
 *
 *  @details
 *  ---------------------------------------------------------------------------------
 *  This header provides:
 *    - A structure (`lim_`) that stores reconstructed variables in each cell.
 *    - Function prototypes for different limiter schemes (e.g., Godunov, Minmod,
 *      Superbee, Weno5, etc.).
 *    - Functions for multi-dimensional reconstructions (1D, 2D, 3D).
 *    - All arrays are now defined using `std::array` with a compile-time constant
 *      (`kEq`) to avoid any form of variable-length arrays.
 *
 *  IMPORTANT NOTES:
 *  ---------------------------------------------------------------------------------
 *  - We assume `eq` is a compile-time macro or constant, e.g., `#define eq 5`, so that
 *    it can be used in a constant expression. We convert it into a named constant
 *    `kEq` (via `enum` or `static constexpr`) to improve clarity and avoid using
 *    "magic numbers" directly in the code.
 *  - If `eq` is not a true compile-time value, you will need to refactor the code
 *    for dynamic allocation (e.g., `std::vector`) instead.
 *  - Using `std::array<double, 2 * kEq>` ensures compile-time checking of
 *    array bounds and conveys the intent more clearly than raw arrays.
 *
 ************************************************************************************/

#ifndef INCLUDE_LIMITERS_HPP_
#define INCLUDE_LIMITERS_HPP_

#include <array>  // for std::array

//----------------------------------------------------------------------------
// We assume 'eq' is a compile-time macro or constant defined elsewhere.
// Here we derive a named compile-time constant 'kEq' from it to use in std::array.
//
// Example:
//   #define eq 5
// or
//   constexpr int eq = 5;
//
// Approach 1 (C++11 or later): static constexpr int kEq = eq;
// Approach 2 (older C++):     enum { kEq = eq };
//
// We'll use the enum approach here for maximum compatibility.
//----------------------------------------------------------------------------
enum { kEq = eq };

//----------------------------------------------------------------------------

/**
 * @brief The structure \b lim_ contains arrays (now `std::array`) for the reconstructed
 *        variables of \f$U\f$ in each cell, stored for each dimension (X1, X2, X3).
 *
 *  - `ux1p`, `ux1m`, etc. store the reconstructed values (plus and minus directions).
 *  - `sx1`, `sx2`, `sx3` might be slopes or auxiliary reconstruction terms.
 *  - `ux` and `U` are additional arrays storing intermediate or final reconstructed data.
 *
 *  The size of each array is defined by the compile-time constant `kEq`:
 *    - `ux1p`, `ux1m`, etc.: each has `2 * kEq` elements.
 *    - `U`: has `kEq + 1` elements, assuming eq+1 is relevant to your physics variables.
 */
struct lim_ {
    std::array<double, 2 * kEq> ux1p;  //!< Reconstructed U in the +X1 direction
    std::array<double, 2 * kEq> ux1m;  //!< Reconstructed U in the -X1 direction
    std::array<double, 2 * kEq> sx1;   //!< Slopes or auxiliary values in the X1 direction

    std::array<double, 2 * kEq> ux2p;  //!< Reconstructed U in the +X2 direction
    std::array<double, 2 * kEq> ux2m;  //!< Reconstructed U in the -X2 direction
    std::array<double, 2 * kEq> sx2;   //!< Slopes or auxiliary values in the X2 direction

    std::array<double, 2 * kEq> ux3p;  //!< Reconstructed U in the +X3 direction
    std::array<double, 2 * kEq> ux3m;  //!< Reconstructed U in the -X3 direction
    std::array<double, 2 * kEq> sx3;   //!< Slopes or auxiliary values in the X3 direction

    std::array<double, 2 * kEq> ux;    //!< A generic or combined set of reconstructed variables
    std::array<double, kEq + 1> U;     //!< Possibly the original or final set of conserved variables
};


/****************************************************************************************
 *                           Function Declarations
 ***************************************************************************************/

/**
 * @brief Driver function for the primitive variable reconstruction.
 *
 * This function likely calls the appropriate reconstruction scheme (e.g., Godunov,
 * Minmod, or WENO) based on the user's choices, and updates the arrays in the `lim_`
 * structure or other global data structures.
 */
void Primitive_Reconstruction();

/**
 * @brief Generic limiter function that dispatches to specific limiter routines
 *        based on the integer `r`.
 *
 * @param A  First argument (usually slope or difference).
 * @param B  Second argument (usually slope or difference).
 * @param r  An integer indicating which limiter to use.
 *
 * @return The limited value.
 */
double Limiter(double A, double B, int r);

/**
 * @brief Godunov limiter function.
 *
 * @param A  First argument.
 * @param B  Second argument.
 *
 * @return The limited value using the Godunov scheme.
 */
double Godunov(double A, double B);

/**
 * @brief Maxmod limiter function.
 *
 * @param A  First argument.
 * @param B  Second argument.
 *
 * @return The limited value using the maxmod scheme.
 */
double Maxmod(double A, double B);

/**
 * @brief Minmod limiter function.
 *
 * @param A  First argument.
 * @param B  Second argument.
 *
 * @return The limited value using the minmod scheme.
 */
double Minmod(double A, double B);

/**
 * @brief MC (monotonized central) limiter function.
 *
 * @param A  First argument.
 * @param B  Second argument.
 *
 * @return The limited value using the MC scheme.
 */
double Mc(double A, double B);

/**
 * @brief Superbee limiter function.
 *
 * @param A  First argument.
 * @param B  Second argument.
 *
 * @return The limited value using the superbee scheme.
 */
double Superbee(double A, double B);

/**
 * @brief WENO5 reconstruction function for 5-point stencils.
 *
 * @param v1  Value at index -2
 * @param v2  Value at index -1
 * @param v3  Value at index  0
 * @param v4  Value at index +1
 * @param v5  Value at index +2
 *
 * @return The reconstructed value using the 5th-order WENO scheme.
 */
double Weno5(double v1, double v2, double v3, double v4, double v5);

/**
 * @brief 1D reconstruction function that takes a pointer to U and a pointer to
 *        a `lim_` structure, plus an index array `I`.
 *
 * @param u  Pointer to the array of variables (likely conserved or primitive).
 * @param l  Pointer to the `lim_` structure where reconstructed values are stored.
 * @param I  Pointer to integer indices specifying the current cell location.
 *
 * @return An integer for success/failure.
 */
int Reconst1D(double *u, lim_ *l, int *I);

/**
 * @brief 2D reconstruction function similar to `Reconst1D` but for 2D.
 *
 * @param u  Pointer to variables.
 * @param l  Pointer to the `lim_` structure.
 * @param I  Pointer to integer indices.
 *
 * @return An integer for success/failure.
 */
int Reconst2D(double *u, lim_ *l, int *I);

/**
 * @brief 3D reconstruction function similar to `Reconst1D` but for 3D.
 *
 * @param u  Pointer to variables.
 * @param l  Pointer to the `lim_` structure.
 * @param I  Pointer to integer indices.
 *
 * @return An integer for success/failure.
 */
int Reconst3D(double *u, lim_ *l, int *I);

#endif  // INCLUDE_LIMITERS_HPP_
