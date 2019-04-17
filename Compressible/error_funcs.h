//==================================================
// Header file containing error functions
//==================================================
#pragma once
#ifndef HEADER_FILE
#define HEADER_FILE
#endif // !HEADER_FILE

#include "print_funcs.h"
#include "constant.h"

vector<vector<double>> error(vector<vector<double>> num, vector<vector<double>> exact); 	// Standard error
double maxChange(vector<vector<double>> prev, vector<vector<double>> curr);	// Maximum change between the previous and current timestep
double L2Norm(vector<vector<double>> error);					// L2 Norm of the error
double L1Norm(vector<vector<double>> error);					// L1 Norm of the error
double Linf(vector<vector<double>> error);					// L infinite Norm of the error


vector<vector<double>> error(vector<vector<double>> num, vector<vector<double>> exact)
{
	vector<vector<double>> err(num.size(), vector<double>(num[0].size()));
	for (int j = 1; j < num.size() - 1; j++)
	{
		for (int i = 1; i < num[0].size() - 1; i++)
		{
			err[j][i] = abs(exact[j][i]) - abs(num[j][i]);
		}
	}
	return err;
}

// Calculates the maximum change of all cells between iterations given the change matrix
double maxChange(vector<vector<double>> prev, vector<vector<double>> curr)
{
	double maxChange = 0;
	vector<vector<double>> delta = error(prev, curr);
	for (int j = 1; j < prev.size() - 1; j++)
	{
		for (int i = 1; i < prev[0].size() - 1; i++)
		{
			delta[j][i] = prev[j][i] - curr[j][i];
			if (abs(delta[j][i]) > maxChange) {
				maxChange = delta[j][i];
			}
		}
	}
	return maxChange;
}

double L2Norm(vector<vector<double>> error)
{
	double L2;
	for (int j = 1; j < error.size() - 1; j++)
	{
		for (int i = 1; i < error[0].size() - 1; i++)
		{
			L2 += pow(error[j][i], 2);
		}
	}
	return sqrt(L2/((error.size() - 2)*(error[0].size() - 2)));
}

double L1Norm(vector<vector<double>> error)
{
	double L1 = 0;
	for (int j = 1; j < error.size() - 1; j++)
	{
		for (int i = 1; i < error[0].size() - 1; i++)
		{
			L1 += abs(error[j][i]);
		}
	}
	return L1 / ((error.size() - 2)*(error[0].size() - 2));
}

double Linf(vector<vector<double>> error)
{
	double Linf = 0;
	for (int j = 1; j < error.size() - 1; j++)
	{
		for (int i = 1; i < error[0].size() - 1; i++)
		{
			if (abs(error[j][i]) > Linf)
			{
				Linf = abs(error[j][i]);
			}
		}
	}
	return Linf;
}
