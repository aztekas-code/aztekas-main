/**
 * @file messages.cpp
 *
 * @brief Helpful functions for aztekas in modern C++ style.
 *
 * This file contains various utility and informational functions:
 *   - Checking system parameters (e.g., number of available cores)
 *   - Managing simulation directories
 *   - Parameter file checks
 *   - Timing and ending messages
 *
 * All C-style I/O has been replaced with C++ stream usage.
 */

#include <iostream>    // For std::cout, std::endl
#include <fstream>     // For std::ofstream
#include <cstdlib>     // For system(), exit
#include <cstdio>      // For popen, pclose
#include <cstring>     // For strcpy, strcat
#include <string>      // For std::string
#include <ctime>       // For clock()
#ifdef _OPENMP
  #include <omp.h>
#endif

#include "main.hpp"    // Includes references to global variables, macros, etc.

//------------------------------------------------------------------------------------
// Check_Sim_Parameters
//------------------------------------------------------------------------------------
void Check_Sim_Parameters() {
#ifdef _OPENMP
    // Obtain the number of physical cores from /proc/cpuinfo (Linux specific)
    FILE* command = popen("grep '^core id' /proc/cpuinfo | sort -u | wc -l", "r");
    if (command) {
        fscanf(command, "%d", &MAX_NUM_THREADS);
        pclose(command);
    }
#endif

    // Using std::cout instead of printf
    std::cout << "\n"
              << "aaaaa  zzzzz  ttttt  eeeee  k   k  aaaaa  sssss\n"
              << "    a     zz    t    e   e  k  k       a  ss   \n"
              << "aaaaa   zzz     t    eeeee  kkk    aaaaa  sssss\n"
              << "a   a  zz       t    e      k  k   a   a     ss\n"
              << "aaaaa  zzzzz    t    eeeee  k   k  aaaaa  sssss\n"
              << "\n"
              << "Running aztekas simulation...\n";

#ifdef _OPENMP
    std::cout << "Parallel version operating with " << OMP_NUM
              << " threads of " << MAX_NUM_THREADS << " available.\n";
#else
    std::cout << "Serial version.\n";
#endif

    std::cout << "\n";

    // Print physics used
    if (PHYSICS == HD)
        std::cout << "Performing a HD simulation.\n";
    if (PHYSICS == RHD)
        std::cout << "Performing a RHD simulation.\n";

    // Equations solved
    std::cout << "Solving the following systems of equations:\n";
#if HYDRO == TRUE
    std::cout << "   Euler inviscid equations (Hyperbolic).\n";
#elif TOV == TRUE
    std::cout << "   Tolman-Oppenheimer-Volkoff (ODE)\n";
#endif

    // Coordinates and metrics
#if PHYSICS == HD
    if (COORDINATES == CARTESIAN)
        std::cout << "Cartesian grid (x,y,z).\n";
    if (COORDINATES == CYLINDRICAL)
        std::cout << "Cylindrical grid (R,z,phi).\n";
  #if POLAR == FALSE
    if (COORDINATES == SPHERICAL)
        std::cout << "Spherical grid (r,theta,phi).\n";
  #elif POLAR == TRUE
    if (COORDINATES == SPHERICAL)
        std::cout << "Polar grid (R,phi).\n";
  #endif

#elif PHYSICS == RHD
  #if COORDINATES == CARTESIAN
    if (METRIC == USER)
        std::cout << "Cartesian grid in a User defined space-time (x,y,z).\n";
    if (METRIC == MINK)
        std::cout << "Cartesian grid in a Minkowski space-time (x,y,z).\n";
    if (METRIC == SCHW)
        std::cout << "Cartesian grid in a Schwarzschild space-time (x,y,z).\n";
    if (METRIC == EF)
        std::cout << "Cartesian grid in an Eddington-Finkelstein space-time (x,y,z).\n";
    if (METRIC == BL)
        std::cout << "Cartesian grid in a Boyer-Lindquist space-time (x,y,z).\n";
    if (METRIC == KS)
        std::cout << "Cartesian grid in a Kerr-Schild space-time (x,y,z).\n";

  #elif COORDINATES == CYLINDRICAL
    if (METRIC == USER)
        std::cout << "Cylindrical grid in a User defined space-time (R,z,phi).\n";
    if (METRIC == MINK)
        std::cout << "Cylindrical grid in a Minkowski space-time (R,z,phi).\n";
    if (METRIC == SCHW)
        std::cout << "Cylindrical grid in a Schwarzschild space-time (R,z,phi).\n";
    if (METRIC == EF)
        std::cout << "Cylindrical grid in an Eddington-Finkelstein space-time (R,z,phi).\n";
    if (METRIC == BL)
        std::cout << "Cylindrical grid in a Boyer-Lindquist space-time (R,z,phi).\n";
    if (METRIC == KS)
        std::cout << "Cylindrical grid in a Kerr-Schild space-time (R,z,phi).\n";

  #elif COORDINATES == SPHERICAL && POLAR == FALSE
    if (METRIC == USER)
        std::cout << "Spherical grid in a User defined space-time (r,theta,phi).\n";
    if (METRIC == MINK)
        std::cout << "Spherical grid in a Minkowski space-time (r,theta,phi).\n";
    if (METRIC == SCHW)
        std::cout << "Spherical grid in a Schwarzschild space-time (r,theta,phi).\n";
    if (METRIC == EF)
        std::cout << "Spherical grid in an Eddington-Finkelstein space-time (r,theta,phi).\n";
    if (METRIC == BL)
        std::cout << "Spherical grid in a Boyer-Lindquist space-time (r,theta,phi).\n";
    if (METRIC == KS)
        std::cout << "Spherical grid in a Kerr-Schild space-time (r,theta,phi).\n";

  #elif COORDINATES == SPHERICAL && POLAR == TRUE
    if (METRIC == USER)
        std::cout << "Polar grid in a User defined space-time (R,phi).\n";
    if (METRIC == MINK)
        std::cout << "Polar grid in a Minkowski space-time (R,phi).\n";
    if (METRIC == SCHW)
        std::cout << "Polar grid in a Schwarzschild space-time (R,phi).\n";
    if (METRIC == EF)
        std::cout << "Polar grid in an Eddington-Finkelstein space-time (R,phi).\n";
    if (METRIC == BL)
        std::cout << "Polar grid in a Boyer-Lindquist space-time (R,phi).\n";
    if (METRIC == KS)
        std::cout << "Polar grid in a Kerr-Schild space-time (R,phi).\n";
  #endif
#endif

    // Resolution
    if (DIM == 1)
        std::cout << "1D simulation with resolution " << Nx1 << " grid cells\n";
    if (DIM == 2)
        std::cout << "2D simulation with resolution " << Nx1 << "X" << Nx2 << " grid cells\n";
    if (DIM == 4)
        std::cout << "2.5D simulation with resolution " << Nx1 << "X" << Nx2 << " grid cells\n";
    if (DIM == 3)
        std::cout << "3D simulation with resolution "
                  << Nx1 << "X" << Nx2 << "X" << Nx3 << " grid cells\n";

    // Equation of state
    if (EOS == IDEAL)
        std::cout << "Ideal equation of state with adiabatic index " << K << ".\n";
    if (EOS == DUST)
        std::cout << "Dust.\n";
    if (EOS == STIFF)
        std::cout << "Stiff equation of state.\n";
    if (EOS == RYU)
        std::cout << "Ryu real relativistic equation of state.\n";

#if HYPERBOLIC == TRUE
    std::cout << "Time integration using a second-order MoL-Runge Kutta.\n";
#endif

#if HYPERBOLIC == TRUE
    // Spatial numerical methods, algorithms, and parameters
    if (RECONST == GODUNOV)
        std::cout << "Zero-order piecewise reconstruction for the primitive variables.\n";
    if (RECONST == MINMOD)
        std::cout << "Second-order piecewise MINMOD reconstruction for the primitive variables.\n";
    if (RECONST == MC)
        std::cout << "Second-order piecewise MC reconstruction for the primitive variables.\n";
    if (RECONST == SUPERBEE)
        std::cout << "Second-order piecewise SUPERBEE reconstruction for the primitive variables.\n";
    if (RECONST == WENO5)
        std::cout << "Fifth-order WENO5 reconstruction for the primitive variables.\n";

    // Flux solver
    if (FLUX == HLL)
        std::cout << "HLL Riemann solver.\n";
    if (FLUX == HLLC)
        std::cout << "HLLC Riemann solver.\n";
#endif

    std::cout << std::endl;

    // ----------------------------------------------------------------
    // Write similar info to the file "info.sim" inside outputdirectory/INFO/
    // ----------------------------------------------------------------
    std::string infoPath = outputdirectory + "INFO/info.sim";
    std::ofstream outFile(infoPath);

    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << infoPath << " for writing.\n";
        return;
    }

    outFile << "\n";
    outFile << "aaaaa  zzzzz  ttttt  eeeee  k   k  aaaaa  sssss\n";
    outFile << "    a     zz    t    e   e  k  k       a  ss   \n";
    outFile << "aaaaa   zzz     t    eeeee  kkk    aaaaa  sssss\n";
    outFile << "a   a  zz       t    e      k  k   a   a     ss\n";
    outFile << "aaaaa  zzzzz    t    eeeee  k   k  aaaaa  sssss\n";
    outFile << "\n";
    outFile << "Running aztekas simulation...\n";

#ifdef _OPENMP
    outFile << "Parallel version operating with " << OMP_NUM
            << " threads of " << MAX_NUM_THREADS << " available.\n";
#else
    outFile << "Serial version.\n";
#endif
    outFile << "\n";

    if (PHYSICS == HD)
        outFile << "Performing a HD simulation.\n";
    if (PHYSICS == RHD)
        outFile << "Performing a RHD simulation.\n";

#if PHYSICS == HD
    if (COORDINATES == CARTESIAN)
        outFile << "Cartesian grid (x,y,z).\n";
    if (COORDINATES == CYLINDRICAL)
        outFile << "Cylindrical grid (R,z,phi).\n";
  #if POLAR == FALSE
    if (COORDINATES == SPHERICAL)
        outFile << "Spherical grid (r,theta,phi).\n";
  #elif POLAR == TRUE
    if (COORDINATES == SPHERICAL)
        outFile << "Polar grid (R,phi).\n";
  #endif

#elif PHYSICS == RHD
  #if COORDINATES == CARTESIAN
    if (METRIC == USER)
        outFile << "Cartesian grid in a User defined space-time (x,y,z).\n";
    if (METRIC == MINK)
        outFile << "Cartesian grid in a Minkowski space-time (x,y,z).\n";
    if (METRIC == SCHW)
        outFile << "Cartesian grid in a Schwarzschild space-time (x,y,z).\n";
    if (METRIC == EF)
        outFile << "Cartesian grid in an Eddington-Finkelstein space-time (x,y,z).\n";
    if (METRIC == BL)
        outFile << "Cartesian grid in a Boyer-Lindquist space-time (x,y,z).\n";
    if (METRIC == KS)
        outFile << "Cartesian grid in a Kerr-Schild space-time (x,y,z).\n";

  #elif COORDINATES == CYLINDRICAL
    if (METRIC == USER)
        outFile << "Cylindrical grid in a User defined space-time (R,z,phi).\n";
    if (METRIC == MINK)
        outFile << "Cylindrical grid in a Minkowski space-time (R,z,phi).\n";
    if (METRIC == SCHW)
        outFile << "Cylindrical grid in a Schwarzschild space-time (R,z,phi).\n";
    if (METRIC == EF)
        outFile << "Cylindrical grid in a Eddington-Finkelstein space-time (R,z,phi).\n";
    if (METRIC == BL)
        outFile << "Cylindrical grid in a Boyer-Lindquist space-time (R,z,phi).\n";
    if (METRIC == KS)
        outFile << "Cylindrical grid in a Kerr-Schild space-time (R,z,phi).\n";

  #elif COORDINATES == SPHERICAL && POLAR == FALSE
    if (METRIC == USER)
        outFile << "Spherical grid in a User defined space-time (r,theta,phi).\n";
    if (METRIC == MINK)
        outFile << "Spherical grid in a Minkowski space-time (r,theta,phi).\n";
    if (METRIC == SCHW)
        outFile << "Spherical grid in a Schwarzschild space-time (r,theta,phi).\n";
    if (METRIC == EF)
        outFile << "Spherical grid in a Eddington-Finkelstein space-time (r,theta,phi).\n";
    if (METRIC == BL)
        outFile << "Spherical grid in a Boyer-Lindquist space-time (r,theta,phi).\n";
    if (METRIC == KS)
        outFile << "Spherical grid in a Kerr-Schild space-time (r,theta,phi).\n";

  #elif COORDINATES == SPHERICAL && POLAR == TRUE
    if (METRIC == USER)
        outFile << "Polar grid in a User defined space-time (R,phi).\n";
    if (METRIC == MINK)
        outFile << "Polar grid in a Minkowski space-time (R,phi).\n";
    if (METRIC == SCHW)
        outFile << "Polar grid in a Schwarzschild space-time (R,phi).\n";
    if (METRIC == EF)
        outFile << "Polar grid in an Eddington-Finkelstein space-time (R,phi).\n";
    if (METRIC == BL)
        outFile << "Polar grid in a Boyer-Lindquist space-time (R,phi).\n";
    if (METRIC == KS)
        outFile << "Polar grid in a Kerr-Schild space-time (R,phi).\n";
  #endif
#endif

    if (DIM == 1)
        outFile << "1D simulation with resolution " << Nx1 << " grid cells\n";
    if (DIM == 2)
        outFile << "2D simulation with resolution " << Nx1 << "X" << Nx2 << " grid cells\n";
    if (DIM == 4)
        outFile << "2.5D simulation with resolution " << Nx1 << "X" << Nx2 << " grid cells\n";
    if (DIM == 3)
        outFile << "3D simulation with resolution "
                << Nx1 << "X" << Nx2 << "X" << Nx3 << " grid cells\n";

    if (EOS == IDEAL)
        outFile << "Ideal equation of state with adiabatic index " << K << ".\n";
    if (EOS == DUST)
        outFile << "Dust.\n";
    if (EOS == STIFF)
        outFile << "Stiff equation of state.\n";
    if (EOS == RYU)
        outFile << "Ryu real relativistic equation of state.\n";

    outFile << "Time integration using a second order MoL-Runge Kutta.\n";

    if (RECONST == GODUNOV)
        outFile << "Using a zero-order piecewise reconstruction for the primitive variables.\n";
    if (RECONST == MINMOD)
        outFile << "Using a second-order piecewise MINMOD reconstruction for the primitive variables.\n";
    if (RECONST == MC)
        outFile << "Using a second-order piecewise MC reconstruction for the primitive variables.\n";
    if (RECONST == SUPERBEE)
        outFile << "Using a second-order piecewise SUPERBEE reconstruction for the primitive variables.\n";
    if (RECONST == WENO5)
        outFile << "Using a fifth-order WENO5 reconstruction for the primitive variables.\n";

    if (FLUX == HLL)
        outFile << "HLL Riemann solver.\n";
    if (FLUX == HLLC)
        outFile << "HLLC Riemann solver.\n";

    outFile.close();
}

//------------------------------------------------------------------------------------
// Manage_Simulation_Info
//------------------------------------------------------------------------------------
void Manage_Simulation_Info(int argc, char* argv[]) {
    // Create outputdirectory
    {
        std::string cmd = "mkdir -p " + outputdirectory;
        system(cmd.c_str());
    }

    // Create the "INFO" subdirectory
    {
        std::string cmd = "mkdir -p " + outputdirectory + "INFO/";
        system(cmd.c_str());
    }

#if MDOT == TRUE
    // Create an "Analysis" subdirectory if needed
    {
        std::string cmd = "mkdir -p " + outputdirectory + "/Analysis/";
        system(cmd.c_str());
    }
#endif

    // Copy paramfile, Makefile, and sources to INFO
    {
        std::string cmd = "cp " + paramfile_name + " " + outputdirectory + "INFO/";
        system(cmd.c_str());
    }
    {
        std::string cmd = "cp Makefile " + outputdirectory + "INFO/";
        system(cmd.c_str());
    }
    {
        std::string cmd = "cp *.c " + outputdirectory + "INFO/";
        system(cmd.c_str());
    }
    {
        std::string cmd = "cp *.h " + outputdirectory + "INFO/";
        system(cmd.c_str());
    }

    Check_Sim_Parameters();

    if (check_param == TRUE)
        std::cin.get();  // Wait for user input (equivalent to getchar())
}

//------------------------------------------------------------------------------------
// Check_Paramfile
//------------------------------------------------------------------------------------
void Check_Paramfile(char* param, int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Wrong number of arguments\n"
                  << "Execute as:\n"
                  << "./aztekas paramfile\n";
        exit(EXIT_FAILURE);
    }

    snprintf(param, argv[1]);
}

//------------------------------------------------------------------------------------
// Computing_Time_Start
//------------------------------------------------------------------------------------
void Computing_Time_Start() {
#ifdef _OPENMP
    start = omp_get_wtime();
    omp_set_num_threads(OMP_NUM);
#else
    start = clock();
#endif
}

//------------------------------------------------------------------------------------
// Ending_Message
//------------------------------------------------------------------------------------
void Ending_Message() {
    std::cout << "\nAZTEKAS termination\n";

#ifdef _OPENMP
    double elapsed = omp_get_wtime() - start;
    if (elapsed > 61.0) {
        int time_sec = static_cast<int>(elapsed);
        int hr  = time_sec / 3600;
        int min = (time_sec % 3600) / 60;
        int sec = (time_sec % 60);
        std::cout << "Expend " << hr << " hr : " << min << " min : "
                  << sec << " sec in the parallelized version using "
                  << OMP_NUM << " threads of " << MAX_NUM_THREADS << " available.\n\n";
    } else {
        std::cout << "Expend " << elapsed << " sec in the parallelized version using "
                  << OMP_NUM << " threads of " << MAX_NUM_THREADS << " available.\n\n";
    }
#else
    double elapsed = static_cast<double>(clock() - start) / CLOCKS_PER_SEC;
    if (elapsed > 61.0) {
        int time_sec = static_cast<int>(elapsed);
        int hr  = time_sec / 3600;
        int min = (time_sec % 3600) / 60;
        int sec = (time_sec % 60);
        std::cout << "Expend " << hr << " hr : " << min << " min : "
                  << sec << " sec in the serial version.\n\n";
    } else {
        std::cout << "Expend " << elapsed << " sec in the serial version.\n\n";
    }
#endif
}
