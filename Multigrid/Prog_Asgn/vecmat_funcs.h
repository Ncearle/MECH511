//=====================================================
// Header file containing matrix manipulation functions
//=====================================================
#pragma once
#ifndef HEADER_FILE
#define HEADER_FILE
#endif // !HEADER_FILE

#include "constant.h"

void AequalB(vector<vector<double> > &A, vector<vector<double> > &B);	// Makes A equal to B
void BequalA(vector<vector<double> > &A, vector<vector<double> > &B);	// Makes A equal to B
void MatAddBtoA(vector<vector<double> > &A, vector<vector<double> > &B);	// A + B w/ no checks
void MatAddAtoB(vector<vector<double> > &A, vector<vector<double> > &B);	// A + B w/ no checks
void resizeMat(vector<vector<double> > &A, int m, int n); // Resize a matrix to m x n
vector<vector<double> > Id(int d);	// Identity matrix
vector<double> MVM(vector<vector<double> > &A, vector<double> &B);	// Matrix vector multiplication
vector<vector<double> > MM(vector<vector<double> > &A, vector<vector<double> > &B);	// Matrix multiplication
vector<double> Vadd(vector<double> &A, vector<double> &B);	// Vector Addition
vector<double> Vsub(vector<double> &A, vector<double> &B);	// Vector Subtraction
vector<vector<double> > Madd(vector<vector<double> > &A, vector<vector<double> > &B);	// Matrix addition
vector<vector<double> > Msub(vector<vector<double> > &A, vector<vector<double> > &B);	// Matrix subtraction
vector<double> ScaV(double S, vector<double> &A);	// Scalar vector multiplication
vector<vector<double> > ScaM(double S, vector<vector<double> > &A);	// Scalar matrix multiplication
vector<vector<vector<double> > > Madd3D(vector<vector<vector<double> > > &A, vector<vector<vector<double> > > &B);	// Matrix addition
vector<vector<vector<double> > > Msub3D(vector<vector<vector<double> > > &A, vector<vector<vector<double> > > &B);	// Matrix subtraction
vector<vector<vector<double> > > ScaM3(double S, vector<vector<vector<double> > > &A); // Scaler 3D matrix multiplication
vector<vector<vector<double> > > copy3(vector<vector<vector<double> > > &A); // Copies a 3D matrix
double MaxV(vector<double> &A);	// Returns the maximum value of a vector

void AequalB(vector<vector<double> > &A, vector<vector<double> > &B)
{
	for (int j = 0; j < A.size(); j++)
	{
		for (int i = 0; i < A[0].size(); i++)
		{
			A[j][i] = B[j][i];
		}
	}
}

void BequalA(vector<vector<double> > &A, vector<vector<double> > &B)
{
	for (int j = 0; j < A.size(); j++)
	{
		for (int i = 0; i < A[0].size(); i++)
		{
			B[j][i] = A[j][i];
		}
	}
}

void MatAddBtoA(vector<vector<double> > &A, vector<vector<double> > &B)
{
	int Br = B.size();
	int Bc = B[0].size();
	for (int r = 0; r < Br; r++)
	{
		for (int c = 0; c < Bc; c++)
		{
			A[r][c] += B[r][c];
		}
	}
}

void MatAddAtoB(vector<vector<double> > &A, vector<vector<double> > &B)
{
	int Br = B.size();
	int Bc = B[0].size();
	for (int r = 0; r < Br; r++)
	{
		for (int c = 0; c < Bc; c++)
		{
			B[r][c] += A[r][c];
		}
	}
}

void resizeMat(vector<vector<double> > &A, int m, int n)
{
	A.resize(m);
	for (int i = 0; i < m; i++)
	{
		A[i].resize(n);
	}
}

vector<vector<double> > Id(int d)
{
	vector<vector<double> > I(d, vector<double>(d));
	for (int j = 0; j < d; j++)
	{
		I[j][j] = 1;
	}
	return I;
}

vector<double> MVM(vector<vector<double> > &A, vector<double> &B)
{
	int Ar = A.size();
	int Ac = A[0].size();
	int Br = B.size();
	vector<double> C(Ar);

	if (Ac == Br)
	{
		for (int r = 0; r < Ar; r++)
		{
			double sum = 0;
			for (int c = 0; c < Br; c++)
			{
				sum += A[r][c] * B[c];
			}
			C[r] = sum;
		}
	}
	else
	{
		cout << "\n!!MULTIPLICATION FAILED -- INNER DIMENSIONS MUST MATCH!!\n";
	}
	return C;
}

vector<vector<double> > MM(vector<vector<double> > &A, vector<vector<double> > &B)
{
	int Ar = A.size();
	int Ac = A[0].size();
	int Br = B.size();
	int Bc = B[0].size();
	vector<vector<double> > C(Ar, vector<double>(Bc));

	if (Ac == Br)
	{
		for (int r = 0; r < Ar; r++)
		{
			for (int c = 0; c < Bc; c++)
			{
				double sum = 0;
				for (int i = 0; i < Ac; i++)
				{
					sum += A[r][i] * B[i][c];
				}
				C[r][c] = sum;
			}
		}
	}
	else
	{
		cout << "\n!!MULTIPLICATION FAILED -- INNER DIMENSIONS MUST MATCH!!\n";
	}

	return C;
}

vector<double> Vadd(vector<double> &A, vector<double> &B)
{
	int Ar = A.size();
	int Br = B.size();
	vector<double> C(Ar);
	if (Ar == Br)
	{
		for (int r = 0; r < Ar; r++)
		{
			C[r] = A[r] + B[r];
		}
	}
	else
	{
		cout << "\n!!ADDITION FAILED -- VECTOR DIMENSIONS MUST MATCH!!\n";
	}
	return C;
}

vector<double> Vsub(vector<double> &A, vector<double> &B)
{
	int Ar = A.size();
	int Br = B.size();
	vector<double> C(Ar);
	if (Ar == Br)
	{
		for (int r = 0; r < Ar; r++)
		{
			C[r] = A[r] - B[r];
		}
	}
	else
	{
		cout << "\n!!SUBTRACTION FAILED -- VECTOR DIMENSIONS MUST MATCH!!\n";
	}
	return C;
}

vector<vector<double> > Madd(vector<vector<double> > &A, vector<vector<double> > &B)
{
	int Ar = A.size();
	int Ac = A[0].size();
	int Bc = B.size();
	int Br = B[0].size();
	vector<vector<double> > C(Ar, vector<double>(Ac));

	if (Ac == Br and Ar == Bc)
	{
		for (int r = 0; r < Ar; r++)
		{
			for (int c = 0; c < Ac; c++)
			{
				C[r][c] = A[r][c] + B[r][c];
			}
		}
	}
	else
	{
		cout << "\n!!ADDITION FAILED -- MATRIX DIMENSIONS MUST MATCH!!\n";
	}
	return C;
}

vector<vector<double> > Msub(vector<vector<double> > &A, vector<vector<double> > &B)
{
	int Ar = A.size();
	int Ac = A[0].size();
	int Bc = B.size();
	int Br = B[0].size();
	vector<vector<double> > C(Ar, vector<double>(Ac));

	if (Ac == Br and Ar == Bc)
	{
		for (int r = 0; r < Ar; r++)
		{
			for (int c = 0; c < Ac; c++)
			{
				C[r][c] = A[r][c] - B[r][c];
			}
		}
	}
	else
	{
		cout << "\n!!SUBTRACTION FAILED -- MATRIX DIMENSIONS MUST MATCH!!\n";
	}
	return C;
}

vector<double> ScaV(double S, vector<double> &A)
{
	int Ar = A.size();
	vector<double> C(Ar);

	for (int r = 0; r < Ar; r++)
	{
		C[r] = S * A[r];
	}
	return C;
}

vector<vector<double> > ScaM(double S, vector<vector<double> > &A)
{
	int Ar = A.size();
	int Ac = A[0].size();

	vector<vector<double> > C(Ar, vector<double>(Ac));
	for (int r = 0; r < Ar; r++)
	{
		for (int c = 0; c < Ac; c++)
		{
			C[r][c] = S * A[r][c];
		}
	}
	return C;
}

vector<vector<vector<double> > > Madd3D(vector<vector<vector<double> > > &A, vector<vector<vector<double> > > &B)
{
	int Ar = A.size();
	int Ac = A[0].size();
	int Ak = A[0][0].size();
	int Bc = B.size();
	int Br = B[0].size();
	int Bk = B[0][0].size();
	vector<vector<vector<double> > > C(Ar, vector<vector<double> >(Ac, vector<double>(Ak)));

	if (Ac == Br and Ar == Bc and Ak == Bk)
	{
		for (int r = 0; r < Ar; r++)
		{
			for (int c = 0; c < Ac; c++)
			{
				for (int k = 0; k < Ak; k++)
				{
					C[r][c][k] = A[r][c][k] + B[r][c][k];
				}
			}
		}
	}
	else
	{
		cout << "\n!!ADDITION FAILED -- MATRIX DIMENSIONS MUST MATCH!!\n";
	}
	return C;
}

vector<vector<vector<double> > > Msub3D(vector<vector<vector<double> > > &A, vector<vector<vector<double> > > &B)
{
	int Ar = A.size();
	int Ac = A[0].size();
	int Ak = A[0][0].size();
	int Bc = B.size();
	int Br = B[0].size();
	int Bk = B[0][0].size();
	vector<vector<vector<double> > > C(Ar, vector<vector<double> >(Ac, vector<double>(Ak)));

	if (Ac == Br and Ar == Bc and Ak == Bk)
	{
		for (int r = 0; r < Ar; r++)
		{
			for (int c = 0; c < Ac; c++)
			{
				for (int k = 0; k < Ak; k++)
				{
					C[r][c][k] = A[r][c][k] - B[r][c][k];
				}
			}
		}
	}
	else
	{
		cout << "\n!!ADDITION FAILED -- MATRIX DIMENSIONS MUST MATCH!!\n";
	}
	return C;
}

vector<vector<vector<double> > > ScaM3(double S, vector<vector<vector<double> > > &A)
{
	int Ar = A.size();
	int Ac = A[0].size();
	int Ak = A[0][0].size();

	vector<vector<vector<double> > > C(Ar, vector<vector<double> >(Ac, vector<double>(Ak)));
	for (int r = 0; r < Ar; r++)
	{
		for (int c = 0; c < Ac; c++)
		{
			for (int k = 0; k < Ak; k++)
			{
				C[r][c][k] = S * A[r][c][k];
			}
		}
	}
	return C;
}

vector<vector<vector<double> > > copy3(vector<vector<vector<double> > > &A)
{
	int Ar = A.size();
	int Ac = A[0].size();
	int Ak = A[0][0].size();

	vector<vector<vector<double> > > C(Ar, vector<vector<double> >(Ac, vector<double>(Ak)));
	for (int r = 0; r < Ar; r++)
	{
		for (int c = 0; c < Ac; c++)
		{
			for (int k = 0; k < Ak; k++)
			{
				C[r][c][k] = A[r][c][k];
			}
		}
	}
	return C;
}

double MaxV(vector<double> &A)
{
	double M = A[0];
	for (int i = 0; i < A.size(); i++)
	{
		if (A[i] > M)
		{
			M = A[i];
		}
	}
	return M;
}
