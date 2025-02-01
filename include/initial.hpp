/************************************************************************************
 *  @file     initial.hpp
 *
 *  @author   Alejandro Aguayo-Ortiz
 *  @brief    Declarations of initial setup functions and variables for aztekas.
 *
 *  @details
 *  ---------------------------------------------------------------------------------
 *  - Declares a global integer `CHECK_NAN` (extern) for tracking NaN checks.
 *  - Provides prototypes for array allocation, mesh resizing, and overall initialization.
 *  - For safer memory management, you may prefer using C++ containers (`std::vector`,
 *    `std::array`) in the corresponding .cpp for `Allocate_Array()`.
 *
 ************************************************************************************/

#ifndef INCLUDE_INITIAL_HPP_
#define INCLUDE_INITIAL_HPP_

/**
 * @brief Global integer used to check for NaN values in the simulation.
 *        Define it in the corresponding .cpp file without 'extern'.
 */
extern int CHECK_NAN;

/**
 * @brief Allocates necessary arrays or data structures for the simulation.
 *        This function could use std::vector or other C++ containers for safety.
 */
void Allocate_Array();

/**
 * @brief Adjusts or allocates mesh/data structures to a new size, if required.
 */
void New_Size();

/**
 * @brief Performs general initialization steps at the start of the simulation.
 */
void Initial();

#endif  // INCLUDE_INITIAL_HPP_
