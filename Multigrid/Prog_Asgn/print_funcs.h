//==================================================
// Header file containing print functions
//==================================================
#pragma once
#ifndef HEADER_FILE
#define HEADER_FILE
#endif // !HEADER_FILE

#include "error_funcs.h"
#include "constant.h"

void printVec(vector<double> vec);		// Prints a vector to the console
void printVec2D(vector<vector<double>> vec);		// Prints a 2D vector to the console
void printVec3D(vector<vector<vector<double>>> vec, int k);		// Prints a layer of a 3D vector to the console 
void vec2D2File(const char* fileName, vector<vector<double>> &vec); // Prints a 2D vector to a file
void vec1D2File(const string fileName, vector<double> &vec); // Prints 1D vector to a file
void printTable(vector<vector<vector<double>>> vec);	// Prints the data in table form 

void printVec(vector<double> vec) 
{
		for (int i = 0; i < vec.size(); i++)
		{
			cout << fixed << setprecision(p) << setw(p+5) << vec[i] << endl;
		}
		cout << endl << endl;
}

void printVec2D(vector<vector<double>> vec) 
{
	for (int j = 0; j < vec.size(); j++)
	{
		for (int i = 0; i < vec[j].size(); i++)
		{
			cout << fixed << setprecision(p) << setw(p+5) << vec[j][i];
		}
		cout << endl;
	}
	cout << endl;
}

void printVec3D(vector<vector<vector<double>>> vec, int k) 
{
	for (int j = vec.size()-1; j > -1; j--)
	{
		for (int i = 0; i < vec[j].size(); i++)
		{
			cout << fixed << setprecision(p) << setw(p+5) << vec[j][i][k];
		}
		cout << endl;
	}
	cout << endl;
}

void vec3D2File(const string fName0, const string fName1, const string fName2, vector<vector<vector<double>>> &vec) // Prints 3D vector to a file
{			
	ofstream outfile0;
	outfile0.open(fName0);
	for (int j = 1; j < vec.size()-1; j++)
	{
		for (int i = 1; i < vec[j].size()-1; i++)
		{
			outfile0 << vec[j][i][0] << " ";		// Separated by a space for .dat files
		}
		outfile0 << endl;
	}
	outfile0 << endl;
	outfile0.close();

	ofstream outfile1;
	outfile1.open(fName1);
	for (int j = 1; j < vec.size()-1; j++)
	{
		for (int i = 1; i < vec[j].size()-1; i++)
		{
			outfile1 << vec[j][i][1] << " ";		// Separated by a space for .dat files
		}
		outfile1 << endl;
	}
	outfile1 << endl;
	outfile1.close();

	ofstream outfile2;
	outfile2.open(fName2);
	for (int j = 1; j < vec.size()-1; j++)
	{
		for (int i = 1; i < vec[j].size()-1; i++)
		{
			outfile2 << vec[j][i][2] << " ";		// Separated by a space for .dat files
		}
		outfile2 << endl;
	}
	outfile2 << endl;
	outfile2.close();
}

void vec2D2File(const string fileName, vector<vector<double>> &vec) // Prints 2D vector to a file
{			
	ofstream outfile;
	outfile.open(fileName);
	for (int j = 0; j < vec.size(); j++)
	{
		for (int i = 0; i < vec[j].size(); i++)
		{
			outfile << vec[j][i] << " ";		// Separated by a space for .dat files
		}
		outfile << endl;
	}
	outfile << endl;
	outfile.close();
}

void vec1D2File(const string fileName, vector<double> &vec) // Prints 1D vector to a file
{
	ofstream outfile;
	outfile.open(fileName);
	for (int i = 0; i < vec.size(); i++)
	{
		outfile << vec[i] << " ";	// Separated by a space for .dat files
	}
	outfile << endl;
	outfile.close();
}

void printTable(vector<vector<vector<double>>> vec)
{
	for (int i = 0; i < imax; i++)
	{
		for (int j = 0; j < jmax; j++)
		{
			cout << "I: " << i;
			cout << "\tJ: " << j;
			cout << "\tP: " << vec[j][i][0];
			cout << "\tu: " << vec[j][i][1];
			cout << "\tv: " << vec[j][i][2];
			cout << endl;
		}
		cout << endl;
	}
}