//Paramfile
int binary;
int check_param;
int restart_simulation, restart_filecount;
char paramfile_name[50], outputdirectory[50], outputfile[50];
char restartfile[50];

int PrintValues(double *tprint, double *dtprint, int *itprint);

int Output1(int *itprint);

int Output2(int *itprint);

int Output3(int *itprint);

int Output1_bin(int *itprint);

int Output2_bin(int *itprint);

int Output3_bin(int *itprint);

void Restart();

void Restart_Bin();

int Read_Parameters_File(char const *paramfile_name);

int User_Parameters(char const *paramfile_name);
