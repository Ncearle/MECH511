//==================================================
// Header file containing exact solution functions
//==================================================
#pragma once
#ifndef HEADER_FILE
#define HEADER_FILE
#endif // !HEADER_FILE

#include "constant.h"

vector<vector<vector<double>>> exactFlux();	// Calculates the exact flux from given values

vector<vector<vector<double>>> exactFlux()
{
	vector<vector<vector<double>>> flux(jmax, vector<vector<double>>(imax, vector<double>(3)));
	for (int j = 0; j < jmax; j++)
	{
		double y = (j - 0.5) / (jmax - 2);
		double Cy = cos(pi*y);
		double Sy = sin(pi*y);
		double C2y = cos(2*pi*y);
		double S2y = sin(2*pi*y);
		for (int i = 0; i < imax; i++)
		{
			double x = (i - 0.5) / (imax - 2);
			double Cx = cos(pi*x);
			double Sx = sin(pi*x);
			double C2x = cos(2*pi*x);
			double S2x = sin(2*pi*x);

			// Pressure
			flux[j][i][0] = -pi/B * (u0*Cx*S2y + v0*S2x*Cy);

			// u velocity
			flux[j][i][1] = P0*pi*Sx*Cy - pow(u0,2)*pi*S2x*pow(S2y,2) - u0*v0*pi*Sx*S2x*(Cy*S2y + 2*C2y*Sy) - u0*5*pow(pi,2)*Sx*S2y/Re;

			// v velocity
			flux[j][i][2] = P0*pi*Cx*Sy - pow(v0,2)*pi*pow(S2x,2)*S2y - u0*v0*pi*Sy*S2y*(Cx*S2x + 2*C2x*Sx) - v0*5*pow(pi,2)*S2x*Sy/Re;
		}
	}
	return flux;
}
