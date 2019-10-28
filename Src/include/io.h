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

void Default_Parameters(char const *paramfile_name);

int User_Parameters(char const *paramfile_name);

void Frequency_Output(double *dtprint);

void Manage_Simulation_Info(int argc, char *argv[]);

void Check_Sim_Parameters();

void Init_Simulation(double *tprint, int *itprint);

<<<<<<< HEAD
int PrintValues(double *tprint, double *dtprint, int *itprint);

int Output1(int *itprint);

int Output2(int *itprint);

int Output3(int *itprint);

int Output1_bin(int *itprint);

int Output2_bin(int *itprint);

int Output3_bin(int *itprint);
=======
void Print_Values(double *tprint, double *dtprint, int *itprint);

void Output1(int *itprint);

void Output2(int *itprint);

void Output3(int *itprint);

void Output1_bin(int *itprint);

void Output2_bin(int *itprint);

void Output3_bin(int *itprint);
>>>>>>> 12b3acd607466560c2bebb7b61677f23252c7907

void Restart();

void Restart_Bin();
