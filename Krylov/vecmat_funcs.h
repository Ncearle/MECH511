//=====================================================
// Header file containing matrix manipulation functions
//=====================================================
#pragma once
#ifndef HEADER_FILE
#define HEADER_FILE
#endif // !HEADER_FILE

#include "constant.h"

void transpose(vector<vector<double> > &A); // Transposes a matrix
void swaprow(vector<vector<double> > &A, int i, int j);	// Swaps row a for b
vector<double> mat2vec(vector<vector<double> > &M);	// Reshapes a matrix into a vector
vector<vector<double> > vec2mat(vector<double> &V);	// Reshapes a vector into a matrix
vector<vector<double> > matrixReshape(vector<vector<double> >& nums, int r, int c); // Reshapes a matrix to another dimension
vector<vector<double> > removeGhost(vector<vector<double> > &M); // Removes ghost cells from a matrix
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
vector<vector<double> > copyMat(vector<vector<double> > &A); // Returns a copy of a matrix
vector<double> copyVec(vector<double> &V); // Returns a copy of a matrix
double MaxV(vector<double> &A);	// Returns the maximum value of a vector

void transpose(vector<vector<double> > &A){

	vector<vector<double> > temp = copyMat(A);
	for(int j = 0; j < A.size(); j++)
	{
		for(int i = 0; i < A[0].size(); i++)
		{
			A[j][i] = temp[i][j];
		}
	}
}
void swaprow(vector<vector<double> > &A, int a, int b){
	vector<double> temp = copyVec(A[a]);
	for (int i = 0; i < A.size(); i++)
	{
			A[a][i] = A[b][i];
			A[b][i] = temp[i];
	}
}
vector<double> mat2vec(vector<vector<double> > &M) {
  int x = M.size();
  int y = M[0].size();

  vector<double> result(x*y);
  int i = 0;
  for (int row = 0; row < x; ++ row)
	{
    for (int col = 0; col < y; ++ col)
		{
      result[i++] = M[row][col];
    }
  }
  return result;
}
vector<vector<double> > vec2mat(vector<double> &V){
	vector<vector<double> > M((jmax-2),vector<double>(imax-2));
	int k = 0;
	for (int j = 0; j < (jmax-2); j++)
	{
		for (int i = 0; i < (imax-2); i++)
		{
			M[j][i] = V[k++];
		}
	}
	return M;
}
vector<vector<double>> matrixReshape(vector<vector<double>>& nums, int r, int c) {
    int x = nums.size();
    int y = nums[0].size();
    if (x * y != r * c) {
        return nums;
    }
    vector<vector<double>> result(r, vector<double>(c));
    int i = 0, j = 0;
    for (int row = 0; row < x; ++ row) {
        for (int col = 0; col < y; ++ col) {
            result[i][j ++] = nums[row][col];
            if (j >= c) {
                i ++;
                j = 0;
            }
        }
    }
    return result;
}
vector<vector<double> > removeGhost(vector<vector<double> > &A){
	int r = A.size()-2;
	int c = A[0].size()-2;
	vector<vector<double> > B((r), vector<double>(c));
	for (int j = 0; j < r; j++)
		for (int i = 0; i < c; i++)
			B[j][i] = A[j+1][i+1];

	return B;
}
void AequalB(vector<vector<double> > &A, vector<vector<double> > &B){
	for (int j = 0; j < A.size(); j++)
	{
		for (int i = 0; i < A[0].size(); i++)
		{
			A[j][i] = B[j][i];
		}
	}
}
void BequalA(vector<vector<double> > &A, vector<vector<double> > &B){
	for (int j = 0; j < A.size(); j++)
	{
		for (int i = 0; i < A[0].size(); i++)
		{
			B[j][i] = A[j][i];
		}
	}
}
void MatAddBtoA(vector<vector<double> > &A, vector<vector<double> > &B){
	int Br = B.size();
	int Bc = B[0].size();
	for (int r = 1; r < Br-1; r++)
	{
		for (int c = 1; c < Bc-1; c++)
		{
			A[r][c] += B[r][c];
		}
	}
}
void MatAddAtoB(vector<vector<double> > &A, vector<vector<double> > &B){
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
void resizeMat(vector<vector<double> > &A, int m, int n){
	A.resize(m);
	for (int i = 0; i < m; i++)
	{
		A[i].resize(n);
	}
}
vector<vector<double> > Id(int d){
	vector<vector<double> > I(d, vector<double>(d));
	for (int j = 0; j < d; j++)
	{
		I[j][j] = 1;
	}
	return I;
}
vector<double> MVM(vector<vector<double> > &A, vector<double> &B){
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
vector<vector<double> > MM(vector<vector<double> > &A, vector<vector<double> > &B){
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
vector<double> Vadd(vector<double> &A, vector<double> &B){
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
vector<double> Vsub(vector<double> &A, vector<double> &B){
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
vector<vector<double> > Madd(vector<vector<double> > &A, vector<vector<double> > &B){
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
vector<vector<double> > Msub(vector<vector<double> > &A, vector<vector<double> > &B){
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
vector<double> ScaV(double S, vector<double> &A){
	int Ar = A.size();
	vector<double> C(Ar);

	for (int r = 0; r < Ar; r++)
	{
		C[r] = S * A[r];
	}
	return C;
}
vector<vector<double> > ScaM(double S, vector<vector<double> > &A){
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
vector<vector<vector<double> > > Madd3D(vector<vector<vector<double> > > &A, vector<vector<vector<double> > > &B){
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
vector<vector<vector<double> > > Msub3D(vector<vector<vector<double> > > &A, vector<vector<vector<double> > > &B){
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
vector<vector<vector<double> > > ScaM3(double S, vector<vector<vector<double> > > &A){
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
vector<vector<vector<double> > > copy3(vector<vector<vector<double> > > &A){
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
vector<vector<double> > copyMat(vector<vector<double> > &A){
	int Ar = A.size();
	int Ac = A[0].size();

	vector<vector<double> > C(Ar, vector<double>(Ac));
	for (int r = 0; r < Ar; r++)
	{
		for (int c = 0; c < Ac; c++)
		{
				C[r][c] = A[r][c];
		}
	}
	return C;
}
vector<double> copyVec(vector<double> &V){
	vector<double> C(V.size());
	for (int r = 0; r < V.size(); r++)
	{
				C[r] = V[r];
	}
	return C;
}
double MaxV(vector<double> &A){
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
