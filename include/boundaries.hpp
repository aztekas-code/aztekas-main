/************************************************************************************
 *  @file     boundaries.hpp
 *
 * .@author   Alejandro Aguayo-Ortiz
 *  @brief    Declarations of boundary condition functions for aztekas.
 *
 *  @details
 *  ---------------------------------------------------------------------------------
 *  This header provides the function prototypes for various boundary condition
 *  implementations (outflow, periodic, reflection, etc.). These functions typically
 *  operate on an array of variables \f$B\f$ representing physical values at the
 *  domain boundaries.
 *
 *  IMPORTANT NOTES:
 *  ---------------------------------------------------------------------------------
 *  - If you want to use modern C++ containers, you could change the function
 *    signatures to accept `std::vector<double>&` or similar.
 *  - Define the implementations of these functions in a .cpp file (e.g., `boundaries.cpp`).
 *
 ************************************************************************************/

#ifndef INCLUDE_BOUNDARIES_HPP_
#define INCLUDE_BOUNDARIES_HPP_

/**
 * @brief Applies the selected boundary conditions to the array \f$B\f$.
 *
 * @param B  Pointer to the array of boundary values.
 */
void Boundaries(double* B);

/**
 * @brief Outflow (or zero-gradient) boundary condition.
 *
 * @param B  Pointer to the array of boundary values.
 */
void Outflow(double* B);

/**
 * @brief Periodic boundary condition, where opposite edges of the domain
 *        match each other's values.
 *
 * @param B  Pointer to the array of boundary values.
 */
void Periodic(double* B);

/**
 * @brief Reflection (or mirror) boundary condition.
 *
 * @param B  Pointer to the array of boundary values.
 */
void Reflection(double* B);

/**
 * @brief User-defined boundary condition, allowing custom logic.
 *
 * @param B  Pointer to the array of boundary values.
 */
void User_Boundaries(double* B);

#endif  // INCLUDE_BOUNDARIES_HPP_
