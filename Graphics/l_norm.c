#include<stdio.h>
#include<math.h>

int main()
{
	FILE *fdata;
	FILE *fexacta;
	char numerical[10];
	char analytic[10];
	int mesh;
   int idum;
   double dum;

   printf("\n");
   printf("################################\n")
   printf("######### L1-NORM - 1D #########");
   printf("################################\n")
   printf("\n");

	printf("Introduce the numerical data\n");
	scanf("%s",&numerical);
	printf("Introduce the analytic date\n");
	scanf("%s",&analytic);

	double Xn[mesh+1];
	double Den_num[mesh+1];
	double Pre_num[mesh+1];
	double Vel_num[mesh+1];

	double Xa[mesh+1];
	double Den_analytic[mesh+1];
	double Pre_analytic[mesh+1];
	double Vel_analytic[mesh+1];

	fdata = fopen(data,"r");
	fexacta = fopen(exacta,"r");

	for(i=0; i<=mesh; i++)
	{
		fscanf(fdata,"%f %f %f %f",&Xd[i],&Nd[i],&Vd[i],&Pd[i]);
		fscanf(fexacta,"%f %f %f %f ",&Xe[i],&Ne[i],&Pe[i],&Ve[i]);
	}

	fclose(fdata);
	fclose(fexacta);

	double Dif[mesh+1];
	double L;

	for(i=0; i<=mesh; i++)
	{
		Dif[i] = fabs(Nd[i] - Ne[i]);
      printf("%d %e \n",i,Dif[i]);
	}
	
	L = 0;

	for(i=0; i<=mesh; i++)
	{
		L = L + Dif[i];
	}

	printf("%f \n",L/mesh);
	return 0;
}
