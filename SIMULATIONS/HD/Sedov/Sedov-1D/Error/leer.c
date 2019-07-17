#include<stdio.h>
#include<math.h>

int main()
{
	FILE *fdata;
	FILE *fexacta;
	int mesh;
	int i;
	char data[10];
	char exacta[10];

	printf("Introduce el numero de puntos:\n");
	scanf("%d",&mesh);
	printf("Introduce el nombre del archivo de datos numericos:\n");
	scanf("%s",&data);
	printf("Introduce el nombre del archivo de datos exactos:\n");
	scanf("%s",&exacta);

	float Xd[mesh+1];
	float Nd[mesh+1];
	float Vd[mesh+1];
	float Pd[mesh+1];

	float Xe[mesh+1];
	float Ne[mesh+1];
	float Ve[mesh+1];
	float Pe[mesh+1];
	float Ee[mesh+1];

	fdata = fopen(data,"r");
	fexacta = fopen(exacta,"r");

	for(i=0; i<=mesh; i++)
	{
		fscanf(fdata,"%f %f %f %f",&Xd[i],&Nd[i],&Vd[i],&Pd[i]);
		fscanf(fexacta,"%f %f %f %f ",&Xe[i],&Ne[i],&Pe[i],&Ve[i]);
	}

	fclose(fdata);
	fclose(fexacta);

	float Dif[mesh+1];
	float L;

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
