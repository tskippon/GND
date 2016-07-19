//

#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>
#include "mex.h"
#include <vector>
#include <cmath>

void eulToBunge(double eul1, double eul2, double eul3, double** bunge);
void transpose(double** a);
void inverse(double** a);
void matrixMultiply(double** a, double** b, double** mult);
double angle(double** m);
void minMisorientation(double** g1, double** g2, double*** symOps, int numSym, double** minMis);
void GNDFunction(int size, double stepx, double stepy, double *phase, double *x, double *y, double *eul1, double *eul2, double *eul3, double *grainId, double *outCurve, double *symType);

extern "C"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//[0] = x, [1] = y, [2] = eul1, [3] = eul2, [4] = eul3, [5] = grainId, [6] = phase, [7] = xstep, [8] = ystep
	mxArray *x_in_m, *y_in_m, *eul1_in_m, *eul2_in_m, *eul3_in_m, *grainId_in_m, *phase_in_m, *symType_in_m;
	mxArray *outCurve_out_m;
	double *x, *y, *eul1, *eul2, *eul3, *phase, *grainId, *symType;
	double* outCurve;
	const mwSize *dims;
	int dimx, dimy, numdims;

	x_in_m = mxDuplicateArray(prhs[0]);
	y_in_m = mxDuplicateArray(prhs[1]);
	eul1_in_m = mxDuplicateArray(prhs[2]);
	eul2_in_m = mxDuplicateArray(prhs[3]);
	eul3_in_m = mxDuplicateArray(prhs[4]);
	grainId_in_m = mxDuplicateArray(prhs[5]);
	phase_in_m = mxDuplicateArray(prhs[6]);
    symType_in_m = mxDuplicateArray(prhs[9]);

	x = mxGetPr(x_in_m);
	y = mxGetPr(y_in_m);
	eul1 = mxGetPr(eul1_in_m);
	eul2 = mxGetPr(eul2_in_m);
	eul3 = mxGetPr(eul3_in_m);
	grainId = mxGetPr(grainId_in_m);
	phase = mxGetPr(phase_in_m);
    symType = mxGetPr(symType_in_m);

	double testval = x[0];
	testval = x[256];

	dims = mxGetDimensions(prhs[0]);
	numdims = mxGetNumberOfDimensions(prhs[0]);
	dimy = (int)dims[0];
	dimx = (int)dims[1];

	double dx = mxGetScalar(prhs[7]);
	double dy = mxGetScalar(prhs[8]);

	outCurve_out_m = plhs[0] = mxCreateDoubleMatrix(dimy, 6, mxREAL);
	outCurve = mxGetPr(outCurve_out_m);
	GNDFunction(dimy, dx, dy, phase, x, y, eul1, eul2, eul3, grainId, outCurve,symType);
}


void eulToBunge(double eul1, double eul2, double eul3, double** bunge)
{
	bunge[0][0] = (std::cos(eul1) * std::cos(eul3)) - (std::sin(eul1) * std::sin(eul3) * std::cos(eul2));
	bunge[1][0] = (std::sin(eul1) * std::cos(eul3)) + (std::cos(eul1) * std::sin(eul3) * std::cos(eul2));
	bunge[2][0] = std::sin(eul3) * std::sin(eul2);
	bunge[0][1] = -(std::cos(eul1) * std::sin(eul3)) - (std::sin(eul1) * std::cos(eul3) * std::cos(eul2));
	bunge[1][1] = -(std::sin(eul1) * std::sin(eul3)) + (std::cos(eul1) * std::cos(eul3) * std::cos(eul2));
	bunge[2][1] = std::cos(eul3) * std::sin(eul2);
	bunge[0][2] = std::sin(eul1) * std::sin(eul2);
	bunge[1][2] = -(std::cos(eul1) * std::sin(eul2));
	bunge[2][2] = std::cos(eul2);
	transpose(bunge);
}

void transpose(double** a)
{
	double holder;
	holder = a[0][1];
	a[0][1] = a[1][0];
	a[1][0] = holder;

	holder = a[0][2];
	a[0][2] = a[2][0];
	a[2][0] = holder;

	holder = a[1][2];
	a[1][2] = a[2][1];
	a[2][1] = holder;
}

void inverse(double** a)
{
	double** inv = new double*[3];
	for (int i = 0; i < 3; ++i)
	{
		inv[i] = new double[3];
	}

	double det = (a[0][0] * a[1][1] * a[2][2]) + (a[1][0] * a[2][1] * a[0][2]) + (a[2][0] * a[0][1] * a[1][2]) - (a[0][0] * a[2][1] * a[1][2]) - (a[1][0] * a[0][1] * a[2][2]) - (a[2][0] * a[1][1] * a[0][2]);

	inv[0][0] = ((a[1][1] * a[2][2]) - (a[1][2] * a[2][1]))/det;
	inv[0][1] = -((a[1][0] * a[2][2]) - (a[1][2] * a[2][0]))/det;
	inv[0][2] = ((a[1][0] * a[2][1]) - (a[1][1] * a[2][0]))/det;
	inv[1][0] = -((a[0][1] * a[2][2]) - (a[0][2] * a[2][1]))/det;
	inv[1][1] = ((a[0][0] * a[2][2]) - (a[0][2] * a[2][0]))/det;
	inv[1][2] = -((a[0][0] * a[2][1]) - (a[0][1] * a[2][0]))/det;
	inv[2][0] = ((a[0][1] * a[1][2]) - (a[0][2] * a[1][1]))/det;
	inv[2][1] = -((a[0][0] * a[1][2]) - (a[0][2] * a[1][0]))/det;
	inv[2][2] = ((a[0][0] * a[1][1]) - (a[0][1] * a[1][0]))/det;

	transpose(inv);

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			a[i][j] = inv[i][j];
		}
	}

	for (int i = 0; i < 3; ++i)
	{
		delete[] inv[i];
	}

	delete[] inv;
}

void matrixMultiply(double** a, double** b, double** mult)
{
	mult[0][0] = a[0][0] * b[0][0] + a[1][0] * b[0][1] + a[2][0] * b[0][2];
	mult[0][1] = a[0][1] * b[0][0] + a[1][1] * b[0][1] + a[2][1] * b[0][2];
	mult[0][2] = a[0][2] * b[0][0] + a[1][2] * b[0][1] + a[2][2] * b[0][2];

	mult[1][0] = a[0][0] * b[1][0] + a[1][0] * b[1][1] + a[2][0] * b[1][2];
	mult[1][1] = a[0][1] * b[1][0] + a[1][1] * b[1][1] + a[2][1] * b[1][2];
	mult[1][2] = a[0][2] * b[1][0] + a[1][2] * b[1][1] + a[2][2] * b[1][2];

	mult[2][0] = a[0][0] * b[2][0] + a[1][0] * b[2][1] + a[2][0] * b[2][2];
	mult[2][1] = a[0][1] * b[2][0] + a[1][1] * b[2][1] + a[2][1] * b[2][2];
	mult[2][2] = a[0][2] * b[2][0] + a[1][2] * b[2][1] + a[2][2] * b[2][2];
}

double angle(double** m)
{
	double ang;
	/*double epsilon;
	double epsilon2;

	epsilon = 0.01;
	epsilon2 = 0.1;

	if (std::abs(m[0][1] - m[1][0]) < epsilon && std::abs(m[0][2] - m[2][0]) < epsilon && std::abs(m[1][2] - m[2][1]) < epsilon)
	{
	if (std::abs(m[0][1] + m[1][0]) < epsilon2 && std::abs(m[0][2] + m[2][0]) < epsilon2 && std::abs(m[1][2] + m[2][1]) < epsilon2 && std::abs(m[0][0] + m[1][1] + m[2][2]) - 3 < epsilon2)
	{
	return 0;
	}
	else
	{
	return 180;
	}
	}*/

	ang = std::acos((m[0][0] + m[1][1] + m[2][2] - 1) / 2) * 180 / M_PI;

	return ang;
}

void minMisorientation(double** g1, double** g2, double*** symOps, int numSym, double** minMis)
{
	transpose(g1);
	matrixMultiply(g1, g2, minMis);
}

//void minMisorientation(double** g1, double** g2, double*** symOps, int numSym, double** minMis)
//{
//	double minMisAng = 1000000;
//	double** mis = new double*[3];
//	double** result = new double*[3];
//	double** tempSym = new double*[3];
//	double ang;
//
//	for (int i = 0; i < 3; ++i)
//	{
//		tempSym[i] = new double[3];
//		result[i] = new double[3];
//		mis[i] = new double[3];
//	}
//
//	for (int n = 0; n < numSym; n++)
//	{
//		for (int i = 0; i < 3; i++)
//		{
//			for (int j = 0; j < 3; j++)
//			{
//				tempSym[i][j] = symOps[n][i][j];
//			}
//		}
//
//		matrixMultiply(g1, tempSym, result);
//		inverse(result);
//		matrixMultiply(result, g2, mis);
//		ang = angle(mis);
//
//		if (ang < minMisAng)
//		{
//			minMisAng = ang;
//			for (int i = 0; i < 3; i++)
//			{
//				for (int j = 0; j < 3; j++)
//				{
//					minMis[i][j] = mis[i][j];
//				}
//			}
//		}
//	}
//
//	for (int i = 0; i < 3; ++i)
//	{
//		delete[] tempSym[i];
//		delete[] mis[i];
//		delete[] result[i];
//	}
//
//	delete[] tempSym;
//	delete[] mis;
//	delete[] result;
//}

void GNDFunction(int size, double stepx, double stepy, double *phase, double *x, double *y, double *eul1, double *eul2, double *eul3, double *grainId, double *outCurve, double *symType) {

	double*** symOpsHex = new double**[12];
	double*** symOpsCube = new double**[24];

	for (int i = 0; i < 12; ++i) {
		symOpsHex[i] = new double*[3];

		for (int j = 0; j < 3; ++j)
			symOpsHex[i][j] = new double[3];
	}

	for (int i = 0; i < 24; ++i) {
		symOpsCube[i] = new double*[3];

		for (int j = 0; j < 3; ++j)
			symOpsCube[i][j] = new double[3];
	}

	double a = 0.5;
	double b = std::sqrt(3) / 2;

	//E //0
	symOpsHex[0][0][0] = 1;
	symOpsHex[0][1][0] = 0;
	symOpsHex[0][2][0] = 0;
	symOpsHex[0][0][1] = 0;
	symOpsHex[0][1][1] = 1;
	symOpsHex[0][2][1] = 0;
	symOpsHex[0][0][2] = 0;
	symOpsHex[0][1][2] = 0;
	symOpsHex[0][2][2] = 1;

	//C6z+ //1
	symOpsHex[1][0][0] = a;
	symOpsHex[1][1][0] = -b;
	symOpsHex[1][2][0] = 0;
	symOpsHex[1][0][1] = b;
	symOpsHex[1][1][1] = a;
	symOpsHex[1][2][1] = 0;
	symOpsHex[1][0][2] = 0;
	symOpsHex[1][1][2] = 0;
	symOpsHex[1][2][2] = 1;

	//C6z- //2
	symOpsHex[2][0][0] = a;
	symOpsHex[2][1][0] = b;
	symOpsHex[2][2][0] = 0;
	symOpsHex[2][0][1] = -b;
	symOpsHex[2][1][1] = a;
	symOpsHex[2][2][1] = 0;
	symOpsHex[2][0][2] = 0;
	symOpsHex[2][1][2] = 0;
	symOpsHex[2][2][2] = 1;

	//C22+ //3
	symOpsHex[3][0][0] = -a;
	symOpsHex[3][1][0] = b;
	symOpsHex[3][2][0] = 0;
	symOpsHex[3][0][1] = b;
	symOpsHex[3][1][1] = a;
	symOpsHex[3][2][1] = 0;
	symOpsHex[3][0][2] = 0;
	symOpsHex[3][1][2] = 0;
	symOpsHex[3][2][2] = -1;

	//C22- //4
	symOpsHex[4][0][0] = a;
	symOpsHex[4][1][0] = -b;
	symOpsHex[4][2][0] = 0;
	symOpsHex[4][0][1] = -b;
	symOpsHex[4][1][1] = -a;
	symOpsHex[4][2][1] = 0;
	symOpsHex[4][0][2] = 0;
	symOpsHex[4][1][2] = 0;
	symOpsHex[4][2][2] = -1;

	//C3z+ //5
	symOpsHex[5][0][0] = -a;
	symOpsHex[5][1][0] = -b;
	symOpsHex[5][2][0] = 0;
	symOpsHex[5][0][1] = b;
	symOpsHex[5][1][1] = -a;
	symOpsHex[5][2][1] = 0;
	symOpsHex[5][0][2] = 0;
	symOpsHex[5][1][2] = 0;
	symOpsHex[5][2][2] = 1;

	//C3z- //6
	symOpsHex[6][0][0] = -a;
	symOpsHex[6][1][0] = b;
	symOpsHex[6][2][0] = 0;
	symOpsHex[6][0][1] = -b;
	symOpsHex[6][1][1] = -a;
	symOpsHex[6][2][1] = 0;
	symOpsHex[6][0][2] = 0;
	symOpsHex[6][1][2] = 0;
	symOpsHex[6][2][2] = 1;

	//C23+ //7
	symOpsHex[7][0][0] = a;
	symOpsHex[7][1][0] = b;
	symOpsHex[7][2][0] = 0;
	symOpsHex[7][0][1] = b;
	symOpsHex[7][1][1] = -a;
	symOpsHex[7][2][1] = 0;
	symOpsHex[7][0][2] = 0;
	symOpsHex[7][1][2] = 0;
	symOpsHex[7][2][2] = -1;

	//C23- //8
	symOpsHex[8][0][0] = -a;
	symOpsHex[8][1][0] = -b;
	symOpsHex[8][2][0] = 0;
	symOpsHex[8][0][1] = -b;
	symOpsHex[8][1][1] = a;
	symOpsHex[8][2][1] = 0;
	symOpsHex[8][0][2] = 0;
	symOpsHex[8][1][2] = 0;
	symOpsHex[8][2][2] = -1;

	//C2z //9
	symOpsHex[9][0][0] = -1;
	symOpsHex[9][1][0] = 0;
	symOpsHex[9][2][0] = 0;
	symOpsHex[9][0][1] = 0;
	symOpsHex[9][1][1] = -1;
	symOpsHex[9][2][1] = 0;
	symOpsHex[9][0][2] = 0;
	symOpsHex[9][1][2] = 0;
	symOpsHex[9][2][2] = 1;

	//C21+ //10
	symOpsHex[10][0][0] = 1;
	symOpsHex[10][1][0] = 0;
	symOpsHex[10][2][0] = 0;
	symOpsHex[10][0][1] = 0;
	symOpsHex[10][1][1] = -1;
	symOpsHex[10][2][1] = 0;
	symOpsHex[10][0][2] = 0;
	symOpsHex[10][1][2] = 0;
	symOpsHex[10][2][2] = -1;

	//C21- //11
	symOpsHex[11][0][0] = -1;
	symOpsHex[11][1][0] = 0;
	symOpsHex[11][2][0] = 0;
	symOpsHex[11][0][1] = 0;
	symOpsHex[11][1][1] = 1;
	symOpsHex[11][2][1] = 0;
	symOpsHex[11][0][2] = 0;
	symOpsHex[11][1][2] = 0;
	symOpsHex[11][2][2] = -1;

	symOpsCube[0][0][0] = 1;
	symOpsCube[0][0][1] = 0;
	symOpsCube[0][0][2] = 0;
	symOpsCube[0][1][0] = 0;
	symOpsCube[0][1][1] = 1;
	symOpsCube[0][1][2] = 0;
	symOpsCube[0][2][0] = 0;
	symOpsCube[0][2][1] = 0;
	symOpsCube[0][2][2] = 1;

	symOpsCube[1][0][0] = 0;
	symOpsCube[1][0][1] = 0;
	symOpsCube[1][0][2] = -1;
	symOpsCube[1][1][0] = 0;
	symOpsCube[1][1][1] = -1;
	symOpsCube[1][1][2] = 0;
	symOpsCube[1][2][0] = -1;
	symOpsCube[1][2][1] = 0;
	symOpsCube[1][2][2] = 0;

	symOpsCube[2][0][0] = 0;
	symOpsCube[2][0][1] = 0;
	symOpsCube[2][0][2] = -1;
	symOpsCube[2][1][0] = 0;
	symOpsCube[2][1][1] = 1;
	symOpsCube[2][1][2] = 0;
	symOpsCube[2][2][0] = 1;
	symOpsCube[2][2][1] = 0;
	symOpsCube[2][2][2] = 0;

	symOpsCube[3][0][0] = -1;
	symOpsCube[3][0][1] = 0;
	symOpsCube[3][0][2] = 0;
	symOpsCube[3][1][0] = 0;
	symOpsCube[3][1][1] = 1;
	symOpsCube[3][1][2] = 0;
	symOpsCube[3][2][0] = 0;
	symOpsCube[3][2][1] = 0;
	symOpsCube[3][2][2] = -1;

	symOpsCube[4][0][0] = 1;
	symOpsCube[4][0][1] = 0;
	symOpsCube[4][0][2] = 0;
	symOpsCube[4][1][0] = 0;
	symOpsCube[4][1][1] = -1;
	symOpsCube[4][1][2] = 0;
	symOpsCube[4][2][0] = 0;
	symOpsCube[4][2][1] = 0;
	symOpsCube[4][2][2] = -1;

	symOpsCube[5][0][0] = -1;
	symOpsCube[5][0][1] = 0;
	symOpsCube[5][0][2] = 0;
	symOpsCube[5][1][0] = 0;
	symOpsCube[5][1][1] = -1;
	symOpsCube[5][1][2] = 0;
	symOpsCube[5][2][0] = 0;
	symOpsCube[5][2][1] = 0;
	symOpsCube[5][2][2] = 1;

	symOpsCube[6][0][0] = 0;
	symOpsCube[6][0][1] = 0;
	symOpsCube[6][0][2] = 1;
	symOpsCube[6][1][0] = 0;
	symOpsCube[6][1][1] = 1;
	symOpsCube[6][1][2] = 0;
	symOpsCube[6][2][0] = -1;
	symOpsCube[6][2][1] = 0;
	symOpsCube[6][2][2] = 0;

	symOpsCube[7][0][0] = 0;
	symOpsCube[7][0][1] = 0;
	symOpsCube[7][0][2] = 1;
	symOpsCube[7][1][0] = 0;
	symOpsCube[7][1][1] = -1;
	symOpsCube[7][1][2] = 0;
	symOpsCube[7][2][0] = 1;
	symOpsCube[7][2][1] = 0;
	symOpsCube[7][2][2] = 0;

	symOpsCube[8][0][0] = 1;
	symOpsCube[8][0][1] = 0;
	symOpsCube[8][0][2] = 0;
	symOpsCube[8][1][0] = 0;
	symOpsCube[8][1][1] = 0;
	symOpsCube[8][1][2] = -1;
	symOpsCube[8][2][0] = 0;
	symOpsCube[8][2][1] = 1;
	symOpsCube[8][2][2] = 0;

	symOpsCube[9][0][0] = 1;
	symOpsCube[9][0][1] = 0;
	symOpsCube[9][0][2] = 0;
	symOpsCube[9][1][0] = 0;
	symOpsCube[9][1][1] = 0;
	symOpsCube[9][1][2] = 1;
	symOpsCube[9][2][0] = 0;
	symOpsCube[9][2][1] = -1;
	symOpsCube[9][2][2] = 0;

	symOpsCube[10][0][0] = -1;
	symOpsCube[10][0][1] = 0;
	symOpsCube[10][0][2] = 0;
	symOpsCube[10][1][0] = 0;
	symOpsCube[10][1][1] = 0;
	symOpsCube[10][1][2] = 1;
	symOpsCube[10][2][0] = 0;
	symOpsCube[10][2][1] = 1;
	symOpsCube[10][2][2] = 0;

	symOpsCube[11][0][0] = -1;
	symOpsCube[11][0][1] = 0;
	symOpsCube[11][0][2] = 0;
	symOpsCube[11][1][0] = 0;
	symOpsCube[11][1][1] = 0;
	symOpsCube[11][1][2] = -1;
	symOpsCube[11][2][0] = 0;
	symOpsCube[11][2][1] = -1;
	symOpsCube[11][2][2] = 0;

	symOpsCube[12][0][0] = 0;
	symOpsCube[12][0][1] = -1;
	symOpsCube[12][0][2] = 0;
	symOpsCube[12][1][0] = 1;
	symOpsCube[12][1][1] = 0;
	symOpsCube[12][1][2] = 0;
	symOpsCube[12][2][0] = 0;
	symOpsCube[12][2][1] = 0;
	symOpsCube[12][2][2] = 1;

	symOpsCube[13][0][0] = 0;
	symOpsCube[13][0][1] = 1;
	symOpsCube[13][0][2] = 0;
	symOpsCube[13][1][0] = -1;
	symOpsCube[13][1][1] = 0;
	symOpsCube[13][1][2] = 0;
	symOpsCube[13][2][0] = 0;
	symOpsCube[13][2][1] = 0;
	symOpsCube[13][2][2] = 1;

	symOpsCube[14][0][0] = 0;
	symOpsCube[14][0][1] = 1;
	symOpsCube[14][0][2] = 0;
	symOpsCube[14][1][0] = 1;
	symOpsCube[14][1][1] = 0;
	symOpsCube[14][1][2] = 0;
	symOpsCube[14][2][0] = 0;
	symOpsCube[14][2][1] = 0;
	symOpsCube[14][2][2] = -1;

	symOpsCube[15][0][0] = 0;
	symOpsCube[15][0][1] = -1;
	symOpsCube[15][0][2] = 0;
	symOpsCube[15][1][0] = -1;
	symOpsCube[15][1][1] = 0;
	symOpsCube[15][1][2] = 0;
	symOpsCube[15][2][0] = 0;
	symOpsCube[15][2][1] = 0;
	symOpsCube[15][2][2] = -1;

	symOpsCube[16][0][0] = 0;
	symOpsCube[16][0][1] = 0;
	symOpsCube[16][0][2] = 1;
	symOpsCube[16][1][0] = 1;
	symOpsCube[16][1][1] = 0;
	symOpsCube[16][1][2] = 0;
	symOpsCube[16][2][0] = 0;
	symOpsCube[16][2][1] = 1;
	symOpsCube[16][2][2] = 0;

	symOpsCube[17][0][0] = 0;
	symOpsCube[17][0][1] = 0;
	symOpsCube[17][0][2] = -1;
	symOpsCube[17][1][0] = -1;
	symOpsCube[17][1][1] = 0;
	symOpsCube[17][1][2] = 0;
	symOpsCube[17][2][0] = 0;
	symOpsCube[17][2][1] = 1;
	symOpsCube[17][2][2] = 0;

	symOpsCube[18][0][0] = 0;
	symOpsCube[18][0][1] = 0;
	symOpsCube[18][0][2] = -1;
	symOpsCube[18][1][0] = 1;
	symOpsCube[18][1][1] = 0;
	symOpsCube[18][1][2] = 0;
	symOpsCube[18][2][0] = 0;
	symOpsCube[18][2][1] = -1;
	symOpsCube[18][2][2] = 0;

	symOpsCube[19][0][0] = 0;
	symOpsCube[19][0][1] = 0;
	symOpsCube[19][0][2] = 1;
	symOpsCube[19][1][0] = -1;
	symOpsCube[19][1][1] = 0;
	symOpsCube[19][1][2] = 0;
	symOpsCube[19][2][0] = 0;
	symOpsCube[19][2][1] = -1;
	symOpsCube[19][2][2] = 0;

	symOpsCube[20][0][0] = 0;
	symOpsCube[20][0][1] = 1;
	symOpsCube[20][0][2] = 0;
	symOpsCube[20][1][0] = 0;
	symOpsCube[20][1][1] = 0;
	symOpsCube[20][1][2] = 1;
	symOpsCube[20][2][0] = 1;
	symOpsCube[20][2][1] = 0;
	symOpsCube[20][2][2] = 0;

	symOpsCube[21][0][0] = 0;
	symOpsCube[21][0][1] = -1;
	symOpsCube[21][0][2] = 0;
	symOpsCube[21][1][0] = 0;
	symOpsCube[21][1][1] = 0;
	symOpsCube[21][1][2] = 1;
	symOpsCube[21][2][0] = -1;
	symOpsCube[21][2][1] = 0;
	symOpsCube[21][2][2] = 0;

	symOpsCube[22][0][0] = 0;
	symOpsCube[22][0][1] = 1;
	symOpsCube[22][0][2] = 0;
	symOpsCube[22][1][0] = 0;
	symOpsCube[22][1][1] = 0;
	symOpsCube[22][1][2] = -1;
	symOpsCube[22][2][0] = -1;
	symOpsCube[22][2][1] = 0;
	symOpsCube[22][2][2] = 0;

	symOpsCube[23][0][0] = 0;
	symOpsCube[23][0][1] = -1;
	symOpsCube[23][0][2] = 0;
	symOpsCube[23][1][0] = 0;
	symOpsCube[23][1][1] = 0;
	symOpsCube[23][1][2] = -1;
	symOpsCube[23][2][0] = 1;
	symOpsCube[23][2][1] = 0;
	symOpsCube[23][2][2] = 0;

	std::vector<double> xoptions;
	std::vector<double> yoptions;

	double minx, maxx;
	double miny, maxy;

	minx = x[0];
	maxx = x[0];
	miny = y[0];
	maxy = y[0];
	xoptions.push_back(x[0]);
	yoptions.push_back(y[0]);

	//check for proper array size
	for (int i = 1; i < size; i++)
	{
		xoptions.push_back(x[i]);
		yoptions.push_back(y[i]);

		if (x[i] < minx)
		{
			minx = x[i];
		}

		if (x[i] > maxx)
		{
			maxx = x[i];
		}

		if (y[i] < miny)
		{
			miny = y[i];
		}

		if (y[i] > maxy)
		{
			maxy = y[i];
		}
	}

	std::sort(xoptions.begin(), xoptions.end());
	std::sort(yoptions.begin(), yoptions.end());
	xoptions.erase(std::unique(xoptions.begin(), xoptions.end()), xoptions.end());
	yoptions.erase(std::unique(yoptions.begin(), yoptions.end()), yoptions.end());

	int dimx = xoptions.size();
	int dimy = yoptions.size();

	//create array
	int** ary = new int*[dimx];
	for (int i = 0; i < dimx; ++i)
	{
		ary[i] = new int[dimy];
	}

	for (int i = 0; i < dimx; i++)
	{
		for (int j = 0; j < dimy; j++)
		{
			ary[i][j] = -1;
		}
	}

	int h, k;
	for (int i = 0; i < size; i++)
	{
		h = std::find(xoptions.begin(), xoptions.end(), x[i]) - xoptions.begin();
		k = std::find(yoptions.begin(), yoptions.end(), y[i]) - yoptions.begin();

		ary[h][k] = i;
	}

	//int h, k;
	////load array with point positions
	//for (int i = 0; i < size; i++)
	//{
	//	h = round((x[i] - minx) / stepx);
	//	k = round((y[i] - miny) / stepy);

	//	ary[h][k] = i;
	//}

	bool nox1 = false;
	bool nox2 = false;
	bool noy1 = false;
	bool noy2 = false;
	double** mis = new double*[3];
	double** g1 = new double*[3];
	double** g2 = new double*[3];
	double rot0;
	double rot1;
	double rot2;
	double dist;

	for (int i = 0; i < 3; i++)
	{
		mis[i] = new double[3];
		g1[i] = new double[3];
		g2[i] = new double[3];
	}

	//compute curves
	for (int i = 0; i < dimx; i++)
	{
		for (int j = 0; j < dimy; j++)
		{
			nox1 = true;
			noy1 = true;
			nox2 = true;
			noy2 = true;

			if (ary[i][j] != -1)
			{
				if (phase[ary[i][j]] != 0)
				{
					int n = 1;
					//for (int n = 1; n <= 2; n++)
					//{
						if (nox1 && i + n < dimx && ary[i + n][j] != -1 && phase[ary[i + n][j]] == phase[ary[i][j]] && grainId[ary[i][j]] == grainId[ary[i + n][j]])
						{
							eulToBunge(eul1[ary[i][j]], eul2[ary[i][j]], eul3[ary[i][j]], g1);
							eulToBunge(eul1[ary[i + n][j]], eul2[ary[i + n][j]], eul3[ary[i + n][j]], g2);

							if (symType[(int)phase[ary[i][j]]-1] == 1)
							{
								minMisorientation(g1, g2, symOpsHex, 12, mis);
							}
							if (symType[(int)phase[ary[i][j]]-1] == 2)
							{
								minMisorientation(g1, g2, symOpsCube, 24, mis);
							}

							rot0 = std::acos((mis[0][0] - 1) / 2);
							rot1 = std::acos((mis[1][1] - 1) / 2);
							rot2 = std::acos((mis[2][2] - 1) / 2);
							dist = std::abs(std::abs(x[ary[i][j]]) - std::abs(x[ary[i + 1][j]]));
							outCurve[0 * size + ary[i][j]] = -((-(mis[1][2]) * rot1 / (2 * std::sin(rot1))) + ((mis[2][1]) * rot2 / (2 * std::sin(rot2)))) / dist;
							outCurve[1 * size + ary[i][j]] = -(((mis[0][2]) * rot0 / (2 * std::sin(rot0))) + (-(mis[2][0]) * rot2 / (2 * std::sin(rot2)))) / dist;
							outCurve[2 * size + ary[i][j]] = -((-(mis[0][1]) * rot0 / (2 * std::sin(rot0))) + ((mis[1][0]) * rot1 / (2 * std::sin(rot1)))) / dist;
							nox1 = false;
						}
					//}

					//for (int n = 1; n <= 2; n++)
					//{
						if (nox2 && nox1 && i - n >= 0 && ary[i - n][j] != -1 && phase[ary[i - n][j]] == phase[ary[i][j]] && grainId[ary[i][j]] == grainId[ary[i - n][j]])
						{
							eulToBunge(eul1[ary[i][j]], eul2[ary[i][j]], eul3[ary[i][j]], g1);
							eulToBunge(eul1[ary[i - n][j]], eul2[ary[i - n][j]], eul3[ary[i - n][j]], g2);

							if (phase[ary[i][j]] == 1)
							{
								minMisorientation(g1, g2, symOpsHex, 12, mis);
							}
							else
							{
								minMisorientation(g1, g2, symOpsCube, 24, mis);
							}

							rot0 = std::acos((mis[0][0] - 1) / 2);
							rot1 = std::acos((mis[1][1] - 1) / 2);
							rot2 = std::acos((mis[2][2] - 1) / 2);
							dist = std::abs(std::abs(x[ary[i][j]]) - std::abs(x[ary[i - n][j]]));
							outCurve[0 * size + ary[i][j]] = -((-(mis[1][2]) * rot1 / (2 * std::sin(rot1))) + ((mis[2][1]) * rot2 / (2 * std::sin(rot2)))) / dist;
							outCurve[1 * size + ary[i][j]] = -(((mis[0][2]) * rot0 / (2 * std::sin(rot0))) + (-(mis[2][0]) * rot2 / (2 * std::sin(rot2)))) / dist;
							outCurve[2 * size + ary[i][j]] = -((-(mis[0][1]) * rot0 / (2 * std::sin(rot0))) + ((mis[1][0]) * rot1 / (2 * std::sin(rot1)))) / dist;
							nox2 = false;
						}
					//}

					//for (int n = 1; n <= 2; n++)
					//{
						if (noy1 && j + n < dimy && ary[i][j + n] != -1 && phase[ary[i][j + n]] == phase[ary[i][j]] && grainId[ary[i][j]] == grainId[ary[i][j + n]])
						{
							eulToBunge(eul1[ary[i][j]], eul2[ary[i][j]], eul3[ary[i][j]], g1);
							eulToBunge(eul1[ary[i][j + n]], eul2[ary[i][j + n]], eul3[ary[i][j + n]], g2);

							if (phase[ary[i][j]] == 1)
							{
								minMisorientation(g1, g2, symOpsHex, 12, mis);
							}
							else
							{
								minMisorientation(g1, g2, symOpsCube, 24, mis);
							}

							rot0 = std::acos((mis[0][0] - 1) / 2);
							rot1 = std::acos((mis[1][1] - 1) / 2);
							rot2 = std::acos((mis[2][2] - 1) / 2);
							dist = std::abs(std::abs(y[ary[i][j]]) - std::abs(y[ary[i][j + n]]));
							outCurve[3 * size + ary[i][j]] = -((-(mis[1][2]) * rot1 / (2 * std::sin(rot1))) + ((mis[2][1]) * rot2 / (2 * std::sin(rot2)))) / dist;
							outCurve[4 * size + ary[i][j]] = -(((mis[0][2]) * rot0 / (2 * std::sin(rot0))) + (-(mis[2][0]) * rot2 / (2 * std::sin(rot2)))) / dist;
							outCurve[5 * size + ary[i][j]] = -((-(mis[0][1]) * rot0 / (2 * std::sin(rot0))) + ((mis[1][0]) * rot1 / (2 * std::sin(rot1)))) / dist;
							noy1 = false;
						}
					//}

					//for (int n = 1; n <= 2; n++)
					//{
						if (noy2 && noy1 && j - n >= 0 && ary[i][j - n] != -1 && phase[ary[i][j - n]] == phase[ary[i][j]] && grainId[ary[i][j]] == grainId[ary[i][j - n]])
						{
							eulToBunge(eul1[ary[i][j]], eul2[ary[i][j]], eul3[ary[i][j]], g1);
							eulToBunge(eul1[ary[i][j - n]], eul2[ary[i][j - n]], eul3[ary[i][j - n]], g2);

							if (phase[ary[i][j]] == 1)
							{
								minMisorientation(g1, g2, symOpsHex, 12, mis);
							}
							else
							{
								minMisorientation(g1, g2, symOpsCube, 24, mis);
							}

							rot0 = std::acos((mis[0][0] - 1) / 2);
							rot1 = std::acos((mis[1][1] - 1) / 2);
							rot2 = std::acos((mis[2][2] - 1) / 2);
							dist = std::abs(std::abs(y[ary[i][j]]) - std::abs(y[ary[i][j - n]]));
							outCurve[3 * size + ary[i][j]] = ((-(mis[1][2]) * rot1 / (2 * std::sin(rot1))) + ((mis[2][1]) * rot2 / (2 * std::sin(rot2)))) / dist;
							outCurve[4 * size + ary[i][j]] = (((mis[0][2]) * rot0 / (2 * std::sin(rot0))) + (-(mis[2][0]) * rot2 / (2 * std::sin(rot2)))) / dist;
							outCurve[5 * size + ary[i][j]] = ((-(mis[0][1]) * rot0 / (2 * std::sin(rot0))) + ((mis[1][0]) * rot1 / (2 * std::sin(rot1)))) / dist;
							noy2 = false;
						}
					//}

					if ((nox1 && nox2) || (noy1 && noy2))
					{
						outCurve[0 * size + ary[i][j]] = 0;
						outCurve[1 * size + ary[i][j]] = 0;
						outCurve[2 * size + ary[i][j]] = 0;
						outCurve[3 * size + ary[i][j]] = 0;
						outCurve[4 * size + ary[i][j]] = 0;
						outCurve[5 * size + ary[i][j]] = 0;
					}
				}
				else
				{
					outCurve[0 * size + ary[i][j]] = 0;
					outCurve[1 * size + ary[i][j]] = 0;
					outCurve[2 * size + ary[i][j]] = 0;
					outCurve[3 * size + ary[i][j]] = 0;
					outCurve[4 * size + ary[i][j]] = 0;
					outCurve[5 * size + ary[i][j]] = 0;
				}
			}
		}
	}

	for (int i = 0; i < 3; ++i)
	{
		delete[] g1[i];
		delete[] mis[i];
		delete[] g2[i];
	}

	delete[] g1;
	delete[] g2;
	delete[] mis;

	for (int i = 0; i < dimx; ++i)
	{
		delete[] ary[i];
	}
	delete[] ary;

	for (int i = 0; i < 12; ++i) {
		for (int j = 0; j < 3; ++j)
			delete[] symOpsHex[i][j];

		delete[] symOpsHex[i];
	}
	delete[] symOpsHex;

	for (int i = 0; i < 24; ++i) {
		for (int j = 0; j < 3; ++j)
			delete[] symOpsCube[i][j];

		delete[] symOpsCube[i];
	}
	delete[] symOpsCube;
}