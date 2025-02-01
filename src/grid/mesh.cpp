/**
 * @file    mesh.cpp
 * @author  Alejandro Aguayo-Ortiz
 *
 * @brief   Cartesian-like mesh grid construction in modern C++.
 *
 * @details
 *  This file replaces the original C version of mesh construction.
 *  It assumes that Nx1, Nx2, Nx3, dx1, dx2, dx3, x1min, x1max, etc. are
 *  defined externally, and that the global `grid` structure uses std::vector<double>
 *  for storing coordinates.
 */

#include "main.hpp"        // or wherever DIM, Nx1, Nx2, Nx3, x1min, etc. are declared
#include "mesh.hpp"        // for grid_ definition, Surface_Volume(), etc.
#include "macros.hpp"      // for any macros like DIM, GRID, HYDRO, etc.

/**
 * @brief Constructs the mesh grid (Cartesian-like).
 *        Sets the dx values and populates grid.X1, grid.X1p, X1m, etc.
 *
 * @return 0 on success, or an error code.
 */
int Mesh() {
    // Local variables
    // M, a, dr, and Delta appear unused except for reference in the original code;
    // keep them if they're needed later or remove if truly unused.
    int i, j, k;
    double M, a, dr, Delta;

#if DIM == 1
    // Compute dx1 using the updated Nx1 (which might already include ghost cells).
    dx1 = (x1max - x1min) / (static_cast<double>(Nx1) - 2.0 * gc);

#elif DIM == 2 || DIM == 4
    dx1 = (x1max - x1min) / (static_cast<double>(Nx1) - 2.0 * gc);
    dx2 = (x2max - x2min) / (static_cast<double>(Nx2) - 2.0 * gc);

#elif DIM == 3
    dx1 = (x1max - x1min) / (static_cast<double>(Nx1) - 2.0 * gc);
    dx2 = (x2max - x2min) / (static_cast<double>(Nx2) - 2.0 * gc);
    dx3 = (x3max - x3min) / (static_cast<double>(Nx3) - 2.0 * gc);
#endif

#if DIM == 1
    // Fill 1D vectors
    for (i = 0; i <= Nx1; i++) {
        grid.X1[i]  = x1min + (i - gc) * dx1;
        grid.X1p[i] = x1min + (i + 0.5 - gc) * dx1;
        grid.X1m[i] = x1min + (i - 0.5 - gc) * dx1;

    #if GRID == LOGMESH
        // Logarithmic mesh example
        grid.X1[i]  = x1min + std::exp(std::log((x1max - x1min + 1.0)) *
                         (i - gc) / (Nx1 - 2.0 * gc)) - 1.0;
        grid.X1p[i] = x1min + std::exp(std::log((x1max - x1min + 1.0)) *
                         (i + 0.5 - gc) / (Nx1 - 2.0 * gc)) - 1.0;
        grid.X1m[i] = x1min + std::exp(std::log((x1max - x1min + 1.0)) *
                         (i - 0.5 - gc) / (Nx1 - 2.0 * gc)) - 1.0;
    #endif
    }

    #if HYDRO == TRUE
      Surface_Volume();
    #endif

#elif DIM == 2 || DIM == 4

    // Fill the X1 arrays
    for (i = 0; i <= Nx1; i++) {
        grid.X1[i]  = x1min + (i - gc) * dx1;
        grid.X1p[i] = x1min + (i + 0.5 - gc) * dx1;
        grid.X1m[i] = x1min + (i - 0.5 - gc) * dx1;

    #if GRID == LOGMESH
        double denom = (Nx1 - 2.0 * gc);
        double factor = std::log((x1max - x1min + lfac) / lfac);
        grid.X1[i]  = x1min + lfac * std::exp(factor * (i - gc) / denom) - lfac;
        grid.X1p[i] = x1min + lfac * std::exp(factor * (i + 0.5 - gc) / denom) - lfac;
        grid.X1m[i] = x1min + lfac * std::exp(factor * (i - 0.5 - gc) / denom) - lfac;
    #endif
    }

    // Fill the X2 arrays
    for (j = 0; j <= Nx2; j++) {
        grid.X2[j]  = x2min + (j - gc) * dx2;
        grid.X2p[j] = x2min + (j + 0.5 - gc) * dx2;
        grid.X2m[j] = x2min + (j - 0.5 - gc) * dx2;
    }

    #if HYDRO == TRUE
      Surface_Volume();
    #endif

#elif DIM == 3

    // Fill the X1 arrays
    for (i = 0; i <= Nx1; i++) {
        grid.X1[i]  = x1min + (i - gc) * dx1;
        grid.X1p[i] = x1min + (i + 0.5 - gc) * dx1;
        grid.X1m[i] = x1min + (i - 0.5 - gc) * dx1;
    }

    // Fill the X2 arrays
    for (j = 0; j <= Nx2; j++) {
        grid.X2[j]  = x2min + (j - gc) * dx2;
        grid.X2p[j] = x2min + (j + 0.5 - gc) * dx2;
        grid.X2m[j] = x2min + (j - 0.5 - gc) * dx2;
    }

    // Fill the X3 arrays
    for (k = 0; k <= Nx3; k++) {
        grid.X3[k]  = x3min + (k - gc) * dx3;
        grid.X3p[k] = x3min + (k + 0.5 - gc) * dx3;
        grid.X3m[k] = x3min + (k - 0.5 - gc) * dx3;
    }

#endif

    return 0;
}
