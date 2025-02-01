/************************************************************************************
 *  @file     main.hpp
 *
 *  @author   Alejandro Aguayo-Ortiz
 *
 *  @brief    Main function header, including necessary C++ standard libraries,
 *            handling OpenMP definitions (if applicable), and declaring global
 *            variables for initialization, integration, boundary conditions,
 *            and I/O routines.
 *
 *  @details
 *  ---------------------------------------------------------------------------------
 *  This header includes:
 *    - C++ standard libraries (iostream, string, vector, cmath, etc.).
 *    - Conditional inclusion for OpenMP if `_OPENMP` is defined (with extern declarations).
 *    - References to internal project headers for mesh, physics, initialization,
 *      integration, boundary conditions, limiters, flux calculations, constants,
 *      macros, I/O, and user parameters.
 *
 *  IMPORTANT NOTES:
 *  ---------------------------------------------------------------------------------
 *  - The code detects if OpenMP is enabled via the `_OPENMP` macro. If so, it
 *    includes `<omp.h>` and defines variables for parallel execution (e.g.,
 *    `MAX_NUM_THREADS`). Otherwise, it uses `<time.h>` for timing.
 *  - This version also uses C++ standard library headers instead of their C
 *    counterparts where possible (`<iostream>` vs. `<stdio.h>`, `<cmath>` vs. `<math.h>`, etc.).
 *  - If you still rely on specific POSIX functions (e.g., `usleep`), `<unistd.h>`
 *    may be necessary. Otherwise, you can remove it.
 *  - The include guard starts with "INCLUDE_AZTEKAS_" and ends with an underscore,
 *    as requested.
 *
 ************************************************************************************/

#ifndef INCLUDE_MAIN_HPP_
#define INCLUDE_MAIN_HPP_

/* Check for OpenMP support */
#ifdef _OPENMP
  #include <omp.h>
  extern int MAX_NUM_THREADS;   //!< Maximum number of threads (OpenMP)
  extern double start;          //!< Start time for performance measurement
#else
  #include <time.h>
  extern clock_t start;         //!< Start time using clock_t if OpenMP is not enabled
#endif

/* Standard C++ libraries */
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstring>

/* If specific POSIX functionality is needed, keep <unistd.h> */
/* #include <unistd.h> */

/* Project headers for mesh and physics (converted to C++ .hpp) */
#include "mesh.hpp"
#include "physics.hpp"

/* Additional project headers for initialization, integration, boundaries, limiters, flux */
#include "initial.hpp"
#include "integration.hpp"
#include "boundaries.hpp"
#include "limiters.hpp"
#include "flux.hpp"

/* Constants, macros, I/O, and user parameters */
#include "macros.hpp"
#include "io.hpp"
#include "user_param.hpp"

#endif // INCLUDE_MAIN_HPP_
