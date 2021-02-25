/**
 * @file io.h
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Input and output function and variable definitions.
 */

// Macros that allow you to define Print_Values0, ..., Print_Values5 for 
// 0, ..., 5 allowed parameters, respectively, and still be called as 
// Print_Values.
#define NARGS(...) NARGS_(0, ## __VA_ARGS__, 5, 4, 3, 2, 1, 0)
#define NARGS_(_5, _4, _3, _2, _1, _0, N, ...) N

#define CONC(A, B) CONC_(A, B)
#define CONC_(A, B) A##B

#define Print_Values(...) CONC(Print_Values_, NARGS(__VA_ARGS__))(__VA_ARGS__)

#define Output_ascii(x) _Generic((x), int: Output_ascii_int, char *: Output_ascii_char)(x);
#define Output_bin(x) _Generic((x), int: Output_bin_int, char *: Output_bin_char)(x);

//Paramfile
int binary;
int numfile;
int check_param;
int restart_simulation, restart_filecount;
double timefile;
char paramfile_name[50], outputdirectory[50], outputfile[50];
char restartfile[50];

void Alternative_Termination();

void Check_Paramfile(char *param, int argc, char *argv[]);

void Check_Sim_Parameters();

void Computing_Time_Start();

void Default_Parameters(char const *paramfile_name);

void Ending_Message();

void Frequency_Output(double *dtprint);

void Init_Simulation(double *tprint, int itprint);

void Manage_Simulation_Info(int argc, char *argv[]);

void Output_ascii_int(int itprint);

void Output_ascii_char(char *itprint);

void Output_bin_int(int itprint);

void Output_bin_char(char *itprint);

void Print_Time_Values(double *tprint, double *dtprint, int itprint);

void Print_Values_0();

void Print_Values_1(char *file_id);

void Restart();

void Restart_Bin();

int User_Parameters(char const *paramfile_name);

void Termination();
