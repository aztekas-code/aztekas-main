#include<stdio.h>
#include<math.h>

int main()
{
	FILE *fdata;
	FILE *fexacta;
	int x1mesh, x2mesh;
	int i,j;
	char data[20];
	char exacta[20];
   double bla;

	printf("Introduce el numero de puntos:\n");
	scanf("%d",&x1mesh);
   x2mesh = 20;
	printf("Introduce el nombre del archivo de datos numericos:\n");
	scanf("%s",&data);
	printf("Introduce el nombre del archivo de datos exactos:\n");
	scanf("%s",&exacta);

	double X1d[(x1mesh+1)*(x2mesh+1)];
	double X2d[(x1mesh+1)*(x2mesh+1)];
	double Nd[(x1mesh+1)*(x2mesh+1)];
	double Pd[(x1mesh+1)*(x2mesh+1)];
	double Ud[(x1mesh+1)*(x2mesh+1)];
	double Vd[(x1mesh+1)*(x2mesh+1)];

	double Xe[(x1mesh+1)];
	double Ne[(x1mesh+1)];
	double Ue[(x1mesh+1)];
	double M[(x1mesh+1)];

	fdata = fopen(data,"r");
	fexacta = fopen(exacta,"r");

	for(i = 0; i <= x1mesh; i++)
	{
      for(j = 0; j <= x2mesh; j++)
      {
         fscanf(fdata,"%lf %lf %lf %lf %lf %lf",&X1d[i*(x2mesh+1) + j],\
         &X2d[i*(x2mesh+1) + j],\
         &Nd[i*(x2mesh+1) + j],\
         &Pd[i*(x2mesh+1) + j],\
         &Ud[i*(x2mesh+1) + j],\
         &Vd[i*(x2mesh+1) + j]);
      }

      fscanf(fexacta,"%lf %lf %lf %lf",&Xe[i],\
      &Ne[i],\
      &Ue[i],\
      &M[i]);
	}

	fclose(fdata);
	fclose(fexacta);

	double Dif[x1mesh];
	double L;

	for(i = 0; i <= x1mesh; i++)
	{
      for(j = x2mesh/2; j <= x2mesh/2; j++)
      {
		   Dif[i] = fabs(Nd[i*(x2mesh + 1) + j] - Ne[x1mesh - i]);
         printf("%e %e %e %e %e \n",X1d[i*(x2mesh+1) + 1],Xe[x1mesh-i],Nd[i*(x2mesh+1) + j],Ne[x1mesh-i],Dif[i]);
      }
	}
	
	L = 0;

	for(i = 0; i <= x1mesh; i++)
	{
		L = L + Dif[i];
	}

   double err = 0;

	for(i=0; i <= x1mesh; i++)
	{
		err = err + fabs(Ne[i]);
	}

	printf("%f \n",L/err);
   
	return 0;
}
