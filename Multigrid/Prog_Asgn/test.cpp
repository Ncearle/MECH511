#include "constant.h"
#include "print_funcs.h"
#include "vecmat_funcs.h"

int main()
{
	vector<vector<double>> A(3, vector<double>(2));
	vector<vector<double>> B(2, vector<double>(3));
	vector<double> Bv(3);
	Bv[0] = 2;
	Bv[1] = 40;
	Bv[2] = 8;


	for (int j = 0; j < A.size(); j++)
	{
		for (int i = 0; i < A[0].size(); i++)
		{
			A[j][i] = (i+1)*(j+1);
		}
	}

	for (int j = 0; j < B.size(); j++)
	{
		for (int i = 0; i < B[0].size(); i++)
		{
			B[j][i] = (i+2)*(j+2);
		}
	}


	// vector<vector<double>> C = Madd(A, B);
	// vector<vector<double>> E = Msub(A, B);
	// vector<double> D = MVM(A, Bv);
	// vector<vector<double>> F = ScaM(1.32445, A);
	vector<vector<double>> G = MM(A,B);

	// printVec(D);

	// printVec2D(C);
	// printVec2D(E);
	// printVec2D(F);
	// printVec2D(G);

	// printVec2D(I(3));

	cout << MaxV(Bv) << endl;

}