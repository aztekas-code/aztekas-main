#include "aztekas/mesh.h"
#include <cmath>

using namespace Mesh;

int Mesh() {
    double dx1 = (x1max - x1min) / (Nx1 - 2 * gc);
    #if DIM >= 2
    double dx2 = (x2max - x2min) / (Nx2 - 2 * gc);
    #endif
    #if DIM == 3
    double dx3 = (x3max - x3min) / (Nx3 - 2 * gc);
    #endif

    for (int i = 0; i <= Nx1; i++) {
        grid.X1[i] = x1min + (i - gc) * dx1;
        grid.X1p[i] = x1min + (i + 0.5 - gc) * dx1;
        grid.X1m[i] = x1min + (i - 0.5 - gc) * dx1;

        #if GRID == LOGMESH
        double logFactor = std::log(x1max - x1min + 1.0);
        grid.X1[i] = x1min + std::exp(logFactor * (i - gc) / (Nx1 - 2 * gc)) - 1;
        grid.X1p[i] = x1min + std::exp(logFactor * (i + 0.5 - gc) / (Nx1 - 2 * gc)) - 1;
        grid.X1m[i] = x1min + std::exp(logFactor * (i - 0.5 - gc) / (Nx1 - 2 * gc)) - 1;
        #endif
    }

    #if DIM >= 2
    for (int j = 0; j <= Nx2; j++) {
        grid.X2[j] = x2min + (j - gc) * dx2;
        grid.X2p[j] = x2min + (j + 0.5 - gc) * dx2;
        grid.X2m[j] = x2min + (j - 0.5 - gc) * dx2;
    }
    #endif

    #if DIM == 3
    for (int k = 0; k <= Nx3; k++) {
        grid.X3[k] = x3min + (k - gc) * dx3;
        grid.X3p[k] = x3min + (k + 0.5 - gc) * dx3;
        grid.X3m[k] = x3min + (k - 0.5 - gc) * dx3;
    }
    #endif

    #if HYDRO == TRUE
    Surface_Volume();
    #endif

    return 0;
}
