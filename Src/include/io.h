/**
 * @file io.h
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Input and output function and variable definitions.
 *
 */

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

void Init_Simulation(double *tprint, int *itprint);

void Manage_Simulation_Info(int argc, char *argv[]);

void Output1(int *itprint);

void Output2(int *itprint);

void Output3(int *itprint);

void Output1_bin(int *itprint);

void Output2_bin(int *itprint);

void Output3_bin(int *itprint);

void Print_Values(double *tprint, double *dtprint, int *itprint);

void Restart();

void Restart_Bin();

int User_Parameters(char const *paramfile_name);

void Termination();
