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
double temperature(vector<double> &Upm);
double pressure (vector<double> &Upm);
vector<double> fluxLimiter(vector<double> &r, string scheme);
vector<double> getPluus(vector<vector<double> > &U, int i, int order, string scheme);
vector<double> getMinus(vector<vector<double> > &U, int i, int order, string scheme);

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
double temperature(vector<double> &Upm){
  double P = pressure(Upm);
  double rho = Upm[0];
  double T = P/(rho*Rair);
  return T;
}
double pressure (vector<double> &Upm){
  double u = uVel(Upm);
  double P = (Upm[2]-0.5*Upm[0]*u*u)*(gam-1);
  return P;
}
vector<double> fluxLimiter(vector<double> &r, string scheme){
  vector<double> FL(3);
  for (int k = 0; k < 3; k++){
    if (scheme == "minmod"){
      vector<double> temp{1, r[k]};
      double x = MinV(temp);
      temp = {0, x};
      FL[k] = MaxV(temp);
    }
    else if (scheme == "vanleer"){
      FL[k] = (r[k] + abs(r[k]))/(1+abs(r[k]));
    }
  }
  return FL;
}

vector<double> getPlus(vector<vector<double> > &U, int i, int order, string scheme=" "){
  vector<double> Uplus(3);
  vector<double> r(3);
  switch (order) {
    case 1: // first order
      for (int k = 0; k < 3; k++){
        Uplus[k] = U[i][k];
      }
    case 2: // second order
      for (int k = 0; k < 3; k++)
      {
        r[k] = (U[i+1][k] - U[i][k])/(U[i][k] - U[i-1][k]);
        if (U[i][k] - U[i-1][k] == 0){r[k] = pow(10, 9);}
      }
      vector<double> FL = fluxLimiter(r, scheme);
      for (int k = 0; k < 3; k++)
      {
        Uplus[k] = U[i][k] + 0.5*FL[k]*(U[i][k] - U[i-1][k]);
      }
  }
  return Uplus;
}
vector<double> getMinus(vector<vector<double> > &U, int i, int order, string scheme=" "){
  vector<double> Uminus(3);
  vector<double> r(3);
  switch (order) {
    case 1: // first order
      for (int k = 0; k < 3; k++){
        Uminus[k] = U[i][k];
      }
    case 2:
      for (int k = 0; k < 3; k++)
      {
        r[k] = (U[i][k] - U[i-1][k])/(U[i+1][k] - U[i][k]);
        if (U[i+1][k] - U[i][k] == 0){r[k] = pow(10, 9);}
      }
      vector<double> FL = fluxLimiter(r, scheme);
      for (int k = 0; k < 3; k++)
      {
        Uminus[k] = U[i][k] + 0.5*FL[k]*(U[i+1][k] - U[i][k]);
      }
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
