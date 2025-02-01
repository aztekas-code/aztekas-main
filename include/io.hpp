/************************************************************************************
 *  @file     io.hpp
 *
 *  @author   Alejandro Aguayo-Ortiz
 *  @brief    Input/Output function and variable definitions for aztekas,
 *adapted to C++.
 *
 *  @details
 *  ---------------------------------------------------------------------------------
 *  - Replaces raw C-style char arrays with std::string for parameters and
 *output file names.
 *  - Retains macros for printing multiple values, preserving the original macro
 *logic.
 *  - Replaces C11 _Generic macros (not valid in standard C++) with inline
 *function overloads.
 *  - All global variables are declared extern; define them in a corresponding
 *.cpp file.
 *
 ************************************************************************************/

#ifndef INCLUDE_IO_HPP_
#define INCLUDE_IO_HPP_

#include <string> // for std::string

/*-----------------------------------------------------------------------------------
 * Macros that allow you to define Print_Values0, ..., Print_Values5 for 0..5
 * parameters, respectively, and still call Print_Values(...) with a varying
 *number of arguments.
 *----------------------------------------------------------------------------------*/
#define NARGS(...) NARGS_(0, ##__VA_ARGS__, 5, 4, 3, 2, 1, 0)
#define NARGS_(_5, _4, _3, _2, _1, _0, N, ...) N
#define CONC(A, B) CONC_(A, B)
#define CONC_(A, B) A##B

#define Print_Values(...) CONC(Print_Values_, NARGS(__VA_ARGS__))(__VA_ARGS__)

/*-----------------------------------------------------------------------------------
 * Output functions: Replaced _Generic with inline function overloads.
 * Simply call Output_ascii(...) or Output_bin(...). Overloads exist for int*
 *and char*.
 *----------------------------------------------------------------------------------*/
inline void
Output_ascii(int *x) { /* Implementation calls Output_ascii_int(x) */
  Output_ascii_int(x);
}
inline void
Output_ascii(char *x) { /* Implementation calls Output_ascii_char(x) */
  Output_ascii_char(x);
}

inline void Output_bin(int *x) { /* Implementation calls Output_bin_int(x) */
  Output_bin_int(x);
}
inline void Output_bin(char *x) { /* Implementation calls Output_bin_char(x) */
  Output_bin_char(x);
}

/*-----------------------------------------------------------------------------------
 * Global variables declared extern.
 * NOTE: Define them (without extern) in your .cpp file (e.g., io.cpp).
 *----------------------------------------------------------------------------------*/
extern int binary;
extern int numfile;
extern int check_param;
extern int restart_simulation;
extern int restart_filecount;
extern double timefile;

/**
 * @brief These file/directory names are now stored as std::string for modern
 * C++ usage. If the code expects C strings elsewhere, use .c_str() when needed.
 */
extern std::string paramfile_name;
extern std::string outputdirectory;
extern std::string outputfile;
extern std::string restartfile;

/*-----------------------------------------------------------------------------------
 * Function Prototypes
 *----------------------------------------------------------------------------------*/
/**
 * @brief Provides an alternative way to terminate the simulation.
 */
void Alternative_Termination();

/**
 * @brief Performs some sort of analysis on array B (likely the main simulation
 * array).
 *
 * @param B pointer to an array of doubles for analysis.
 */
void Analysis(double *B);

/**
 * @brief Reads and checks parameter file arguments from the command line.
 *
 * @param param  Pointer to the parameter name or buffer.
 * @param argc   Argument count from main().
 * @param argv   Argument vector from main().
 */
void Check_Paramfile(char *param, int argc, char *argv[]);

/**
 * @brief Checks internal simulation parameters after reading config.
 */
void Check_Sim_Parameters();

/**
 * @brief Starts a timer or records a timestamp for computing time measurement.
 */
void Computing_Time_Start();

/**
 * @brief Loads default parameters from a given parameter file name.
 *
 * @param paramfile_name  Name of the parameter file.
 */
void Default_Parameters(char const *paramfile_name);

/**
 * @brief Prints or logs a final message before the program ends.
 */
void Ending_Message();

/**
 * @brief Manages output frequency based on dtprint.
 *
 * @param dtprint Pointer to the time interval for output.
 */
void Frequency_Output(double *dtprint);

/**
 * @brief Initializes the simulation, possibly reading from config or setting
 * initial data.
 *
 * @param tprint  Pointer to the current output time.
 * @param itprint Pointer to the current output iteration.
 */
void Init_Simulation(double *tprint, int *itprint);

/**
 * @brief Handles simulation info, likely reading command line args, etc.
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 */
void Manage_Simulation_Info(int argc, char *argv[]);

/**
 * @brief Actual (legacy) function that writes ASCII output given an int
 * pointer.
 */
void Output_ascii_int(int *itprint);

/**
 * @brief Actual (legacy) function that writes ASCII output given a char
 * pointer.
 */
void Output_ascii_char(char *itprint);

/**
 * @brief Actual (legacy) function that writes binary output given an int
 * pointer.
 */
void Output_bin_int(int *itprint);

/**
 * @brief Actual (legacy) function that writes binary output given a char
 * pointer.
 */
void Output_bin_char(char *itprint);

/**
 * @brief Prints time-related values (e.g., current simulation time) to screen
 * or file.
 *
 * @param tprint   Pointer to the current output time.
 * @param dtprint  Pointer to the time interval for output.
 * @param itprint  Pointer to the current output iteration.
 */
void Print_Time_Values(double *tprint, double *dtprint, int *itprint);

/**
 * @brief Handles printing with no parameters.
 *        Called by the macro Print_Values_0().
 */
void Print_Values_0();

/**
 * @brief Handles printing with one parameter (char*).
 *        Called by the macro Print_Values_1().
 *
 * @param file_id  A string or file identifier.
 */
void Print_Values_1(char *file_id);

/**
 * @brief Restarts the simulation from a stored state.
 */
void Restart();

/**
 * @brief Restarts the simulation from a binary-stored state.
 */
void Restart_Bin();

/**
 * @brief Reads and sets user-defined parameters from the specified file.
 *
 * @param paramfile_name  Name of the parameter file.
 * @return An integer code indicating success/failure.
 */
int User_Parameters(char const *paramfile_name);

/**
 * @brief Terminates the simulation, possibly freeing resources or printing
 * final status.
 */
void Termination();

#endif // INCLUDE_IO_HPP_
