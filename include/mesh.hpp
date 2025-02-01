/************************************************************************************
 *  @file     mesh.hpp
 *
 *  @author   Alejandro Aguayo-Ortiz
 *  @brief    Definitions of the vectors, structures, and functions needed for
 *            the mesh grid at the center and at the interfaces.
 *
 *  @details
 *  ---------------------------------------------------------------------------------
 *  This header provides:
 *    - Global variables for the number of grid points (Nx1, Nx2, Nx3) and the
 *      domain boundaries (x1min, x1max, etc.).
 *    - A C++ structure `grid_` (adapted with `std::vector`) to hold the mesh
 *      information, including positions at cell centers and interfaces,
 *      and additional surfaces (S1p, S1m, etc.).
 *    - C++ structures `gauge_` and `der_gauge_` that store metric components
 *      and their derivatives, now using `std::array` for better safety.
 *    - Prototypes for functions related to mesh construction, surface volume
 *      calculation, and time-stepping.
 *
 *  IMPORTANT NOTES:
 *  ---------------------------------------------------------------------------------
 *  - All original comments have been translated into English, retaining and
 *    expanding them when necessary.
 *  - Global variables remain `extern` to preserve linkage across translation units.
 *  - Where possible, C-style arrays have been replaced by `std::vector` and `std::array`
 *    to leverage C++ safety and expressiveness.
 *
 ************************************************************************************/

#ifndef INCLUDE_MESH_HPP_
#define INCLUDE_MESH_HPP_

#include <vector>
#include <array>

/**
 * @brief Defines the number of grids in each dimension.
 *        These are declared as extern for linkage across modules.
 */
extern int Nx1, Nx2, Nx3;

/**
 * @brief Defines the domain boundaries in each dimension (min and max values).
 *        These are declared as extern for linkage across modules.
 */
extern double x1max, x2max, x3max;
extern double x1min, x2min, x3min;

/**
 * @brief Structure representing the mesh grid and its properties.
 *
 *  - `time`: current simulation time.
 *  - `X1`, `X2`, `X3`: vectors for cell-center coordinates in each dimension.
 *  - `X1p`, `X1m`, etc.: vectors for positive/negative interface coordinates.
 *  - `S1p`, `S1m`, etc.: vectors for surface or area elements on interfaces.
 *
 *  This structure is adapted to C++ by using `std::vector<double>` instead of
 *  raw pointers. Memory management is thus handled automatically, and the size
 *  of each vector can be set according to Nx1, Nx2, Nx3 at runtime.
 */
struct grid_ {
    double time;              //!< Current simulation time

    std::vector<double> X1;   //!< Cell-center coordinates in the X1 dimension
    std::vector<double> X1p;  //!< Positive interface coordinates in the X1 dimension
    std::vector<double> X1m;  //!< Negative interface coordinates in the X1 dimension

    std::vector<double> X2;   //!< Cell-center coordinates in the X2 dimension
    std::vector<double> X2p;  //!< Positive interface coordinates in the X2 dimension
    std::vector<double> X2m;  //!< Negative interface coordinates in the X2 dimension

    std::vector<double> X3;   //!< Cell-center coordinates in the X3 dimension
    std::vector<double> X3p;  //!< Positive interface coordinates in the X3 dimension
    std::vector<double> X3m;  //!< Negative interface coordinates in the X3 dimension

    std::vector<double> S1p;  //!< Surface (or area) for the positive side of X1
    std::vector<double> S1m;  //!< Surface (or area) for the negative side of X1

    std::vector<double> S2p;  //!< Surface (or area) for the positive side of X2
    std::vector<double> S2m;  //!< Surface (or area) for the negative side of X2

    std::vector<double> S3p;  //!< Surface (or area) for the positive side of X3
    std::vector<double> S3m;  //!< Surface (or area) for the negative side of X3
};

/**
 * @brief Global instance of the mesh grid.
 *
 * You can resize the vectors inside `grid` according to Nx1, Nx2, Nx3
 * during the initialization phase in `Mesh()`.
 */
extern grid_ grid;

/**
 * @brief Structure storing metric information (lapse, shift, and metric tensors)
 *        at a particular location in space, for general relativity or other curved
 *        space-time formulations.
 *
 *  - `I[3]`: integer indices in the 3D mesh.
 *  - `x[4]`: coordinates (e.g., time + 3 spatial coords).
 *  - `lapse`: the lapse function.
 *  - `beta_con` / `beta_cov`: contravariant/covariant shift vectors.
 *  - `gamma_con` / `gamma_cov`: 3D metric (contravariant/covariant).
 *  - `g_con` / `g_cov`: 4D metric (contravariant/covariant).
 *  - `dety`: determinant (or related factor) of the spatial metric.
 */
struct gauge_ {
    std::array<int, 3>    I;         //!< (i, j, k) mesh indices
    std::array<double, 4> x;         //!< (t, x, y, z) or similar 4D coordinates

    double lapse;                    //!< Lapse function

    std::array<double, 3> beta_con;  //!< Contravariant shift vector
    std::array<double, 3> beta_cov;  //!< Covariant shift vector

    std::array<std::array<double, 3>, 3> gamma_con;  //!< 3D contravariant metric
    std::array<std::array<double, 3>, 3> gamma_cov;  //!< 3D covariant metric

    std::array<std::array<double, 4>, 4> g_con;     //!< 4D contravariant metric
    std::array<std::array<double, 4>, 4> g_cov;     //!< 4D covariant metric

    double dety;                    //!< Determinant or related factor of the metric
};

/**
 * @brief Structure storing derivatives of metric components (or gauge functions).
 *
 *  - `dlapse`: partial derivatives of the lapse function in the 3 spatial directions.
 *  - `dbeta`: partial derivatives of the shift vector in the 3 spatial directions.
 *  - `dgam`: partial derivatives of the 3D metric components.
 */
struct der_gauge_ {
    std::array<double, 3> dlapse;  //!< dlapse/dx, dlapse/dy, dlapse/dz
    std::array<std::array<double, 3>, 3> dbeta;     //!< derivative of beta vector
    std::array<std::array<std::array<double, 3>, 3>, 3> dgam;  //!< derivative of gamma metric
};


/**
 * @brief Spatial resolutions (grid spacing) and related parameters.
 *        Declared as extern for usage throughout the program.
 */
extern double dx1, dx2, dx3;  //!< Spatial step sizes in each dimension
extern double dt;             //!< Time step
extern double tmax;           //!< Maximum time
extern double cou;            //!< CFL or Courant number for stability

/**
 * @brief Creates or initializes the mesh.
 *        This function likely allocates arrays or resizes vectors and sets domain values.
 *
 * @return An integer indicating success (0) or failure.
 */
int Mesh();

/**
 * @brief Calculates surface or volume elements (areas in 2D or 3D) for the interfaces.
 *        The function name suggests it computes surface volumes (perhaps face areas
 *        or cell volumes).
 */
void Surface_Volume();

/**
 * @brief Computes the time step (dt) based on mesh resolutions and the Courant condition.
 *
 * @return The computed time step (dt).
 */
double TimeStep();

#endif  // INCLUDE_MESH_HPP_
