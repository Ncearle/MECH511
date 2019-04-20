//==========================================================
// Header file containing physical functions and other stuff
//==========================================================
#pragma once
#ifndef HEADER_FILE
#define HEADER_FILE
#endif // !HEADER_FILE

#include "constant.h"

double speedofsound(vector<double> &Upm);
double uVel(vector<double> &Upm);
double density(vector<double> &Upm);
// double temperature(vector<double> &Upm);
double pressure (vector<double> &Upm);
vector<double> getPlus(vector<vector<double> > &U, int i);
vector<double> getMinus(vector<vector<double> > &U, int i);
vector<vector<double> > getF(vector<vector<double> > &U);

double speedofsound(vector<double> &Upm){
  double P = pressure(Upm);
  double c = sqrt(gam*P/Upm[0]);
  return c;
}
double uVel(vector<double> &Upm){
  return Upm[1]/Upm[0];
}
double density(vector<double> &Upm){
  return Upm[0];
}
// double temperature(vector<double> &Upm){
//   double u = Upm[1]/Upm[0];
//   double rho = Upm[0];
//   double T = (Upm[2] - rho*u*u/2)/(rho*Cv);
//   return T;
// }
double pressure (vector<double> &Upm){
  double u = uVel(Upm);
  double P = (Upm[2]-0.5*Upm[0]*u*u)*(gam-1);
  return P;
}

vector<double> getPlus(vector<vector<double> > &U, int i){
  vector<double> Uplus(3);
  for (int k = 0; k < 3; k++)
  {
    Uplus[k] = U[i][k]; // first order
    // Uplus[k] = (3*U[i][k] - U[i-1][k])/2; //second order
  }
  return Uplus;
}
vector<double> getMinus(vector<vector<double> > &U, int i){
  vector<double> Uminus(3);
  for (int k = 0; k < 3; k++)
  {
    Uminus[k] = U[i][k]; // first order
    // Uminus[k] = (3*U[i][k] - U[i+1][k])/2; //second order
  }
  return Uminus;
}

vector<vector<double> > getF(vector<vector<double> > &U){
  vector<vector<double> > F(imax, vector<double>(3));
  for (int i = 0; i < imax; i++)
  {
    double u = uVel(U[i]);
    double P = pressure(U[i]);

    F[i][0] = U[i][0] * u;
    F[i][1] = U[i][1] * u + P;
    F[i][2] = U[i][2] * u + P*u;
  }
  return F;
}
