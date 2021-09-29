#include <stdio.h> 		             /*Import Standard Module*/
#include <stdlib.h> 		           /*Import Standard Library*/
#include <math.h>		               /*Import Math Module*/
#include <complex.h>          	   /*Import Modules for Complex Numbers*/
#include "header.h"		             /*Import Module for Dawson Function*/


complex double Phi[501][501][81];
complex double Phi2[501][101][81];
complex double S[501][501][81];
complex double K[501][501][81];
complex double rho[501][81];

double PNorm[81];
double Prob[501][81];


int main() {

complex double P;
double R;


double Kpot;
double lambda;
double wnull;


Kpot=1.0;
lambda=0.1;
wnull=1;



double h;
int tmax;
int lmax;
int itr, itr1, itr2, itr3, itr4;



h=0.1;
lmax=10/h+1;
tmax=8/h+1;

double x[lmax],xf[lmax];
double t[tmax];


for (itr = 0; itr < lmax; itr++) {
x[itr]= -5+itr*h;
xf[itr]= -5+itr*h;
}

for (itr = 0; itr < tmax; itr++) {
t[itr]= itr*h;
}




/*Phase Integrals*/


for (itr1 = 0; itr1 < lmax; itr1++) {
  for (itr2 = 0; itr2 < lmax; itr2++) {
    for (itr3 = 0; itr3 < tmax; itr3++) {

      for (itr4 =8 ; itr4 < itr3; itr4++) {
        P= P+h*omega(x[itr1], x[itr2], t[itr4]);
      }
      Phi[itr1][itr2][itr3]=P;
      P=0;
    }
  }
}


for (itr1 = 0; itr1 < lmax; itr1++) {
  for (itr2 = 0; itr2 < lmax; itr2++) {
    for (itr3 = 0; itr3 < tmax; itr3++) {

      for (itr4 =8 ; itr4 < itr3; itr4++) {
        P= P+h/(omega(x[itr1], x[itr2], t[itr4])*omega(x[itr1], x[itr2], t[itr4]));
      }
      Phi2[itr1][itr2][itr3]=P;
      P=0;
    }
  }
}


/*Path Integral*/

for (itr1 = 0; itr1 < lmax; itr1++) {
  for (itr2 = 0; itr2 < lmax; itr2++) {
    for (itr3 = 0; itr3 < tmax; itr3++) {

      K[itr1][itr2][itr3]=csqrt(1/(2*M_PI*I*t[itr3]))*cexp(I*Scl(x[itr1], x[itr2], t[itr3]))*(1+(I*lambda*Kpot)*Phi2[itr1][itr2][itr3]*csqrt(Phi[itr1][itr2][itr3]/csin(Phi[itr1][itr2][itr3]))) ;
    }
  }
}



/*Integrating the wave function*/


for (itr1 = 0; itr1 < lmax; itr1++) {
    for (itr3 = 0; itr3 < tmax; itr3++) {
      for (itr2 = 0; itr2 < lmax; itr2++) {
        P= P+h*K[itr1][itr2][itr3]*Wpkt(x[itr2]);
      }
      rho[itr1][itr3]=P;
      P=0;
  }
}

/*Computing Proability*/


R=0;

for (itr3 = 9; itr3 < tmax; itr3++) {
    for (itr1 = 0; itr1 < 50; itr1++) {
      R=R+cabsf(rho[itr1][itr3])*cabsf(rho[itr1][itr3]);   
  }
    for (itr1 = 51; itr1 < lmax; itr1++) {
      R=R+cabsf(rho[itr1][itr3])*cabsf(rho[itr1][itr3]);   
  }
  PNorm[itr3]=R;
  R=0;
}




for (itr3 = 9; itr3 < tmax; itr3++) {
    for (itr1 = 0; itr1 < lmax; itr1++) {
      Prob[itr1][itr3]=cabsf(rho[itr1][itr3])*cabsf(rho[itr1][itr3])/PNorm[itr3];
    }
}


/*Plot Statements*/

/*
for (itr1 = 0; itr1 < lmax; itr1++) {
  for (itr2 = 0; itr2 < tmax; itr2++) {
        printf("%f\t%f\t%f\n", x[itr1], t[itr2],  crealf(csqrt(Phi[itr1][20][itr2]/csin(Phi[itr1][20][itr2]))));
  }
}


for (itr1 = 0; itr1 < lmax; itr1++) {
  for (itr2 = 0; itr2 < tmax; itr2++) {
        printf("%f\t%f\t%f\n", x[itr1], t[itr2],  cimagf((Phi2[itr1][20][itr2])));
  }
}
*/
for (itr1 = 0; itr1 < lmax; itr1++) {
  for (itr2 = 9; itr2 < tmax; itr2++) {
        printf("%f\t%f\t%f\n", x[itr1],  t[itr2],  Prob[itr1][itr2] );
  }
}
/*
for (itr2 = 9; itr2 < tmax; itr2++) {
      printf("%f\t%f\n",  t[itr2],  PNorm[itr2] );
}
*/

}
