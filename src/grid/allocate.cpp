/**
 * @file allocate.cpp
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Essential allocation functions for aztekas using modern C++ containers.
 *
 * @details
 *  This file replaces the original C-style memory allocations (malloc/free) with
 *  C++-style `std::vector<double>` resizing. The global arrays and grid vectors
 *  are assumed to be declared as `std::vector<double>` in the appropriate headers.
 */

#include <vector>
#include "mesh.hpp"         // For grid_ structure and Nx1, Nx2, Nx3
#include "integration.hpp"  // For U, Q, U0, U1, etc.
#include "initial.hpp"      // For Allocate_Array(), New_Size()
#include "macros.hpp"       // For DIM, eq, etc. if needed
#include "user_param.h"     // Potentially for gc (ghost cells) or other parameters

/**
 * @brief Adjusts the global Nx1, Nx2, Nx3 to account for ghost cells.
 *        E.g., Nx1 = Nx1 + 2*gc, etc.
 *
 * This version is unchanged from the original, just placed in a .cpp file.
 */
void New_Size() {
    Nx1 = Nx1 + 2 * gc;
    Nx2 = Nx2 + 2 * gc;
    Nx3 = Nx3 + 2 * gc;
}

/**
 * @brief Allocates/resizes all vectors used in aztekas according to DIM and eq.
 *
 * This modern C++ version uses `std::vector<double>::resize()` instead of `malloc()`.
 * The total sizes depend on Nx1, Nx2, Nx3 (which are updated by `New_Size()`),
 * and the number of equations eq (or kEq).
 */
void Allocate_Array() {
    // First, resize Nx1, Nx2, Nx3 to include ghost cells
    New_Size();

#if DIM == 1
    // Resize the 1D grid vectors
    grid.X1 .resize(Nx1 + 1);
    grid.X1p.resize(Nx1 + 1);
    grid.X1m.resize(Nx1 + 1);

    grid.S1p.resize(Nx1 + 1);
    grid.S1m.resize(Nx1 + 1);

    // Main solution arrays (size = (Nx1+1)*(eq+1))
    U.resize((Nx1 + 1) * (eq + 1));
    Q.resize((Nx1 + 1) * (eq + 1));

    // Time-stepping auxiliaries
    U0.resize((Nx1 + 1) * (eq + 1));
    U1.resize((Nx1 + 1) * (eq + 1));
    U2.resize((Nx1 + 1) * (eq + 1));
    Q0.resize((Nx1 + 1) * (eq + 1));
    Q1.resize((Nx1 + 1) * (eq + 1));
    Q2.resize((Nx1 + 1) * (eq + 1));

    // Reconstructed interface variables
    U1p.resize((Nx1 + 1) * (eq + 1));
    U1m.resize((Nx1 + 1) * (eq + 1));

#elif DIM == 2 || DIM == 4
    // Resize the 2D grid vectors
    grid.X1 .resize(Nx1 + 1);
    grid.X1p.resize(Nx1 + 1);
    grid.X1m.resize(Nx1 + 1);
    grid.X2 .resize(Nx2 + 1);
    grid.X2p.resize(Nx2 + 1);
    grid.X2m.resize(Nx2 + 1);

    // Interfaces
    grid.S1p.resize((Nx1 + 1) * (Nx2 + 1));
    grid.S1m.resize((Nx1 + 1) * (Nx2 + 1));
    grid.S2p.resize((Nx1 + 1) * (Nx2 + 1));
    grid.S2m.resize((Nx1 + 1) * (Nx2 + 1));

    // Main solution arrays (size = (Nx1+1)*(Nx2+1)*(eq+1))
    U.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));
    Q.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));

    // Time-stepping auxiliaries
    U0.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));
    U1.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));
    U2.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));
    Q0.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));
    Q1.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));
    Q2.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));

    // Reconstructed interface variables
    U1p.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));
    U1m.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));
    U2p.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));
    U2m.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));

#elif DIM == 3
    // Resize the 3D grid vectors
    grid.X1 .resize(Nx1 + 1);
    grid.X1p.resize(Nx1 + 1);
    grid.X1m.resize(Nx1 + 1);
    grid.X2 .resize(Nx2 + 1);
    grid.X2p.resize(Nx2 + 1);
    grid.X2m.resize(Nx2 + 1);
    grid.X3 .resize(Nx3 + 1);
    grid.X3p.resize(Nx3 + 1);
    grid.X3m.resize(Nx3 + 1);

    // Interfaces
    grid.S1p.resize((Nx1 + 1) * (Nx2 + 1) * (Nx3 + 1));
    grid.S1m.resize((Nx1 + 1) * (Nx2 + 1) * (Nx3 + 1));
    grid.S2p.resize((Nx1 + 1) * (Nx2 + 1) * (Nx3 + 1));
    grid.S2m.resize((Nx1 + 1) * (Nx2 + 1) * (Nx3 + 1));
    grid.S3p.resize((Nx1 + 1) * (Nx2 + 1) * (Nx3 + 1));
    grid.S3m.resize((Nx1 + 1) * (Nx2 + 1) * (Nx3 + 1));

    // Main solution array
    U.resize((Nx1 + 1) * (Nx2 + 1) * (Nx3 + 1) * (eq + 1));

    // Time-stepping auxiliaries for U, Q, etc.
    // (Note: the original code did not allocate Q in 3D, but you may need it.)
    U0.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));
    U1.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));
    U2.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));
    Q0.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));
    Q1.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));
    Q2.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));

    // Reconstructed interface variables
    U1p.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));
    U1m.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));
    U2p.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));
    U2m.resize((Nx1 + 1) * (Nx2 + 1) * (eq + 1));
#endif
}
