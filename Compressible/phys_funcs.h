//=====================================================
// Header file containing physical functions
//=====================================================
#pragma once
#ifndef HEADER_FILE
#define HEADER_FILE
#endif // !HEADER_FILE

#include "constant.h"

double inittemp(double rho, double P);
double speedofsound(vector<vector<double> > &U, int i)
double uVel(vector<vector<dobule> > &U, int i);
vector<double> Uplus(vector<vector<double> > &U, int i);
vector<double> Uminus(vector<vector<double> > &U, int i);

double inittemp(double rho, double P){
  double T = P/(rho*R);
  return T;
}
double speedofsound(vector<vector<double> > &U, int i){
  double T = (U[2][i] - pow(U[1][i],2)/(2*U[0][i])) / (U[0][i]*Cv);
  double c = sqrt(gam*R*T);
  return c;
}
double uVel(vector<vector<double> > &U, int i, string PM){
  double u;
  if (PM == "+"){
    u = U[1][i]/U[0][i];    // first order u_i+1/2;i = u_i
  }
  else if (PM == "-"){
    u = U[1][i]/U[0][i];
  }
  return u;
}
double density(vector<vector<double> > &U, int i, string PM){
  double rho
  if (PM == "+"){
    rho = U[0][i];   // first order rho_i+1/2;i = rho_i
  }
  else if (PM == "-"){
    rho = U[0][i];
  }
  return rho;
}

vector<double> Uplus(vector<vector<double> > &U, int i){
  vector<double> Uplus(3);
  for (int k = 0; k < 3; k++)
  {
    Uplus[k] = U[k][i]; // first order
    // Uplus[k] = (3*U[k][i] - U[k][i-1])/2; //second order
  }
  return Uplus;
}
vector<double> Uminus(vector<vector<double> > &U, int i){
  vector<double> Uminus(3);
  for (int k = 0; k < 3; k++)
  {
    Uminus[k] = U[k][i]; // first order
    // Uminurs[k] = (3*U[k][i] - U[k][i-1])/2; //second order
  }
  return Uminus;
}
