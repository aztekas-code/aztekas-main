/**
 * @file flux.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Numerical flux computing and implementation.
 */

//Do not erase any of these libraries//
#include"main.h"

int Flux1D(vec_ *v, lim_ *l, int *I)
{
   flx_ f;

   AMATRIX1D(l->ux,v,I);
#if riemann == 1 //Hll
   VECTOR(1,'f',l,&f,I);
   Hll(v->Fp,&f,1);
   VECTOR(0,'f',l,&f,I);
   Hll(v->Fm,&f,1);
#elif riemann == 2//Hllc
   VECTOR(1,'f',l,&f,I);
   Hllc(v->Fp,&f,1);
   VECTOR(0,'f',l,&f,I);
   Hllc(v->Fm,&f,1);
#endif

   return 0;
}

int Flux2D(vec_ *v, lim_ *l, int *I)
{
   flx_ f;

   AMATRIX2D(l->ux,v,I);
#if riemann == 1//Hll
   VECTOR(1,'f',l,&f,I);
   Hll(v->Fp,&f,1);
   VECTOR(1,'g',l,&f,I);
   Hll(v->Gp,&f,2);
   VECTOR(0,'f',l,&f,I);
   Hll(v->Fm,&f,1);
   VECTOR(0,'g',l,&f,I);
   Hll(v->Gm,&f,2);
#elif riemann == 2//Hllc
   VECTOR(1,'f',l,&f,I);
   Hllc(v->Fp,&f,1);
   VECTOR(1,'g',l,&f,I);
   Hllc(v->Gp,&f,2);
   VECTOR(0,'f',l,&f,I);
   Hllc(v->Fm,&f,1);
   VECTOR(0,'g',l,&f,I);
   Hllc(v->Gm,&f,2);
#endif

   return 0;
}

int Flux3D(vec_ *v, lim_ *l, int *I)
{
   int n;

   flx_ f;

   AMATRIX3D(l->ux,v,I);
   VECTOR(1,'f',l,&f,I);
   Hll(v->Fp,&f,1);
   VECTOR(1,'g',l,&f,I);
   Hll(v->Gp,&f,2);
   VECTOR(1,'h',l,&f,I);
   Hll(v->Hp,&f,3);
   VECTOR(0,'f',l,&f,I);
   Hll(v->Fm,&f,1);
   VECTOR(0,'g',l,&f,I);
   Hll(v->Gm,&f,2);
   VECTOR(0,'h',l,&f,I);
   Hll(v->Hm,&f,3);

   return 0;
}

int Hll(double *F, flx_ *f, int x)
{
   int n;
   double q;
   double sR, sL;
   double QR[eq+1], QL[eq+1];
   double FR[eq+1], FL[eq+1];
   double qp[eq+1], qm[eq+1];
   double fp[eq+1], fm[eq+1];

   sR = f->lp;
   sL = f->lm;

   if(sL >= 0)
   {
      for(n = 0; n < eq; n++)
      {
         FL[n] = f->fm[n];
         F[n] = FL[n];
      }
   }
   else if(sL <= 0 && sR >= 0)
   {
      for(n = 0; n < eq; n++)
      {
         QR[n] = f->qp[n];
         QL[n] = f->qm[n];
         FR[n] = f->fp[n];
         FL[n] = f->fm[n];
         q = sR*sL*(QR[n] - QL[n]);
         F[n] = (sR*FL[n] - sL*FR[n] + q)/(sR - sL);
      }
   }
   else if(f->lp <= 0)
   {
      for(n = 0; n < eq; n++)
      {
         FR[n] = f->fp[n];
         F[n] = FR[n];
      }
   }

   return 0;
}

int Hllc(double *F, flx_ *f, int x)
{
   int n;
   double pstar, ustar, rhobar, abar;
   double rhoR, rhoL, pR, pL;
   double uR, uL, vR, vL, wR, wL;
   double aR, aL;
   double sR, s, sL;
   double qR, qL, eR, eL;
   double QR[eq+1], QL[eq+1];
   double FR[eq+1], FL[eq+1];
   double QsR[eq+1], QsL[eq+1];
   double D[eq+1];
   double rsu;

   rhoR = f->up[0];
   rhoL = f->um[0];
   pR   = f->up[1];
   pL   = f->um[1];
#if DIM == 1
   uR   = f->up[2];
   uL   = f->um[2];
   vR   = 0.0;
   vL   = 0.0;
   wR   = 0.0;
   wL   = 0.0;
#elif DIM == 2
   uR   = f->up[2];
   uL   = f->um[2];
   vR   = f->up[3];
   vL   = f->um[3];
   wR   = 0.0;
   wL   = 0.0;
#elif DIM == 3 || DIM == 4
   uR   = f->up[2];
   uL   = f->um[2];
   vR   = f->up[3];
   vL   = f->um[3];
   wR   = f->up[4];
   wL   = f->um[4];
#endif

   if(x == 1)
   {
      aR = sqrt(K*pR/rhoR);
      aL = sqrt(K*pL/rhoL);

      rhobar = 0.5*(rhoR + rhoL);
      abar   = 0.5*(aR + aL);

      pstar = 0.5*(pR + pL) + 0.5*(uL - uR)*(rhobar*abar);

      if(pstar <= pL)
      {
         qL = 1.0;
      }
      else
      {
         qL = sqrt(1.0 + ((K + 1.0)/(2.0*K))*((pstar/pL) - 1.0));
      }

      if(pstar <= pR)
      {
         qR = 1.0;
      }
      else
      {
         qR = sqrt(1.0 + ((K + 1.0)/(2.0*K))*((pstar/pR) - 1.0));
      }

      sR = uR + aR*qR;
      sL = uL - aL*qL;
      ustar = (pR - pL + rhoL*uL*(sL - uL) - rhoR*uR*(sR - uR))/(rhoL*(sL - uL) - rhoR*(sR - uR));
      s = ustar;
      
      eR = 0.5*rhoR*(uR*uR + vR*vR + wR*wR) + pR/(K-1.0);
      eL = 0.5*rhoL*(uL*uL + vL*vL + wL*wL) + pL/(K-1.0);

      rsu = rhoR*((sR - uR)/(sR - s));

      QsR[0] = rsu;
      QsR[1] = rsu*(eR/rhoR + (s - uR)*(s + pR/(rhoR*(sR - uR))));
      QsR[2] = rsu*s;
      QsR[3] = rsu*vR;
      QsR[4] = rsu*wR;

      rsu = rhoL*((sL - uL)/(sL - s));

      QsL[0] = rsu;
      QsL[1] = rsu*(eL/rhoL + (s - uL)*(s + pL/(rhoL*(sL - uL))));
      QsL[2] = rsu*s;
      QsL[3] = rsu*vL;
      QsL[4] = rsu*wL;
   }

   if(x == 2)
   {
      aR = sqrt(K*pR/rhoR);
      aL = sqrt(K*pL/rhoL);

      rhobar = 0.5*(rhoR + rhoL);
      abar   = 0.5*(aR + aL);

      pstar = 0.5*(pR + pL) + 0.5*(vL - vR)*(rhobar*abar);

      if(pstar <= pL)
      {
         qL = 1.0;
      }
      else
      {
         qL = sqrt(1 + ((K + 1)/(2*K))*((pstar/pL) - 1.0));
      }

      if(pstar <= pR)
      {
         qR = 1.0;
      }
      else
      {
         qR = sqrt(1 + ((K + 1)/(2*K))*((pstar/pR) - 1.0));
      }

      sR = vR + aR*qR;
      sL = vL - aL*qL;
      ustar = (pR - pL + rhoL*vL*(sL - vL) - rhoR*vR*(sR - vR))/(rhoL*(sL - vL) - rhoR*(sR - vR));
      s = ustar;
      
      eR = 0.5*rhoR*(uR*uR + vR*vR + wR*wR) + pR/(K-1.0);
      eL = 0.5*rhoL*(uL*uL + vL*vL + wL*wL) + pL/(K-1.0);

      rsu = rhoR*((sR - vR)/(sR - s));

      QsR[0] = rsu;
      QsR[1] = rsu*(eR/rhoR + (s - vR)*(s + pR/(rhoR*(sR - vR))));
      QsR[2] = rsu*uR;
      QsR[3] = rsu*s;
      QsR[4] = rsu*wR;

      rsu = rhoL*((sL - vL)/(sL - s));

      QsL[0] = rsu;
      QsL[1] = rsu*(eL/rhoL + (s - vL)*(s + pL/(rhoL*(sL - vL))));
      QsL[2] = rsu*uL;
      QsL[3] = rsu*s;
      QsL[4] = rsu*wL;
   }

   if(sL >= 0)
   {
      for(n = 0; n < eq; n++)
      {
         F[n] = f->fm[n];
      }
   }
   else if(sL < 0 && s >= 0)
   {
      for(n = 0; n < eq; n++)
      {
         FL[n] = f->fm[n];
         QL[n] = f->qm[n];
         F[n] = FL[n] + sL*(QsL[n] - QL[n]);
      }
   }
   else if(s <= 0 && sR >= 0)
   {
      for(n = 0; n < eq; n++)
      {
         FR[n] = f->fp[n];
         QR[n] = f->qp[n];
         F[n] = FR[n] + sR*(QsR[n] - QR[n]);
      }
   }
   else if(sR <= 0)
   {
      for(n = 0; n < eq; n++)
      {
         F[n] = f->fp[n];
      }
   }

   return 0;
}
