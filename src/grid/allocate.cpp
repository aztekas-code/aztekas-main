#include "aztekas/mesh.h"

using namespace Mesh;

void New_Size() {
    Nx1 += 2 * gc;
    Nx2 += 2 * gc;
    Nx3 += 2 * gc;
}

void Allocate_Array(int eq, int gc) {
    New_Size();

    int size1D = (Nx1 + 1);
    int size2D = (Nx1 + 1) * (Nx2 + 1);
    int size3D = (Nx1 + 1) * (Nx2 + 1) * (Nx3 + 1);

    grid = Grid(size1D, size2D, size3D);

    int solutionSize = (Nx1 + 1) * (eq + 1);
    #if DIM >= 2
    solutionSize *= (Nx2 + 1);
    #endif
    #if DIM == 3
    solutionSize *= (Nx3 + 1);
    #endif

    U.resize(solutionSize);
    Q.resize(solutionSize);
    U0.resize(solutionSize);
    U1.resize(solutionSize);
    U2.resize(solutionSize);
    Q0.resize(solutionSize);
    Q1.resize(solutionSize);
    Q2.resize(solutionSize);
    U1p.resize(solutionSize);
    U1m.resize(solutionSize);

    #if DIM >= 2
    U2p.resize(solutionSize);
    U2m.resize(solutionSize);
    #endif

    #if DIM == 3
    U3p.resize(solutionSize);
    U3m.resize(solutionSize);
    #endif
}
