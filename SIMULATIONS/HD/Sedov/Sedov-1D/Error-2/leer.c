#include<stdio.h>
#include<math.h>

int main()
{
	FILE *fdata;
	FILE *fexacta;
	int x1mesh, x2mesh;
	int i,j;
	char data[10];
	char exacta[10];
   double bla;

	printf("Introduce el numero de puntos:\n");
	scanf("%d",&x1mesh);
	printf("Introduce el numero de puntos:\n");
	scanf("%d",&x2mesh);
	printf("Introduce el nombre del archivo de datos numericos:\n");
	scanf("%s",&data);
	printf("Introduce el nombre del archivo de datos exactos:\n");
	scanf("%s",&exacta);

	double X1d[x1mesh+1];
	double X2d[x2mesh+1];
	double Nd[x1mesh*x2mesh];
	double Vd[x1mesh*x2mesh];
	double Pd[x1mesh*x2mesh];

	double Xe[x1mesh+1];
	double Ne[x1mesh+1];
	double Ve[x1mesh+1];
	double Pe[x1mesh+1];
	double Ee[x1mesh+1];

	fdata = fopen(data,"r");
	fexacta = fopen(exacta,"r");

	for(i=0; i<=x1mesh; i++)
	{
	   for(j=0; j<=x2mesh; j++)
	   {
		   fscanf(fdata,"%lf %lf %lf %lf %lf %lf",&X1d[i],&X2d[j],&Nd[i*x2mesh + j],&bla,&bla,&bla);
      }
		fscanf(fexacta,"%lf %lf %lf %lf ",&Xe[i],&Ne[i],&Pe[i],&Ve[i]);
	}

	fclose(fdata);
	fclose(fexacta);

	double Dif[x1mesh+1];
	double L;

	for(i=0; i<=x1mesh; i++)
	{
		Dif[i] = fabs(Nd[i*x2mesh + x2mesh/2] - Ne[i]);
      printf("%d %e %e %e \n",i,Nd[i*x2mesh + x2mesh/2],Ne[i],Dif[i]);
	}
	
	L = 0;

	for(i=0; i<=x1mesh; i++)
	{
		L = L + Dif[i];
	}

	printf("%f \n",L/x1mesh);
	return 0;
}
