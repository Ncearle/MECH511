//=====================================================
// Header file containing physical functions
//=====================================================
#pragma once
#ifndef HEADER_FILE
#define HEADER_FILE
#endif // !HEADER_FILE

#include "constant.h"

double inittemp(double rho, double P);
double speedofsound(vector<double> &Upm);
double uVel(vector<double> &Upm);
double density(vector<double> &Upm);
vector<double> Uplus(vector<vector<double> > &U, int i);
vector<double> Uminus(vector<vector<double> > &U, int i);

double inittemp(double rho, double P){
  double T = P/(rho*R);
  return T;
}
double speedofsound(vector<double> &Upm){
  double T = (Upm[2] - pow(Upm[1],2)/(2*Upm[0])) / (Upm[0]*Cv);
  double c = sqrt(gam*R*T);
  return c;
}
double uVel(vector<double> &Upm){
  return Upm[1]/Upm[0];
}
double density(vector<double> &Upm){
  return Upm[0];
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
    // Uminus[k] = (3*U[k][i] - U[k][i+1])/2; //second order
  }
  return Uminus;
}
