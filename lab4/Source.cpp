//lab4, v20. Simon Karasik, course 2, group 5.
#include <fstream>
#include <string>
#include <numeric>
#include "MatrixMath.h"

using namespace std;

const double EPS = 10E-5;

struct MethodSummary
{
	Matrix B;
	Vector b;
	Vector error;
	Vector x;
	int stepCnt;
	int theoreticalStepCnt;
};

bool eq(double x, double y) {
	return abs(x - y) <= EPS;
}

bool eq(const Vector & x, const Vector & y)
{
	if (x.size() != y.size())
		throw invalid_argument("Vectors have different size");
	for (int i = 0; i < x.size(); i++)
		if (!eq(x[i], y[i]))
			return false;
	return true;
}

double computeNorm(const Vector & v)
{
	double norm = 0;
	for (int i = 0; i < v.size(); i++)
		norm = max(norm, abs(v[i]));
	return norm;
}

double computeNorm(const Matrix & mat) 
{
	double norm = 0;
	for (int i = 0; i < getN(mat); i++)
	{
		double sum = 0;
		for (int j = 0; j < getM(mat); j++)
			sum += abs(mat[i][j]);
		norm = max(norm, sum);
	}

	return norm;
}

double computeCovergenceSpeed(const Matrix & B, const Vector & b)
{
	double normB = computeNorm(B), normb = computeNorm(b);
	double numerator = log(EPS*(1 - normB) / normb);
	double denumerator = log(normB);
	return trunc(numerator / denumerator) + 1;
}


void printMethodSummary(
	ostream & os, 
	const string & methodName, 
	const MethodSummary & methodSummary)
{
	os << methodName << endl;
	if (!methodSummary.B.empty())
	{ 
		os << "B:" << endl;
		printMatrix(os, methodSummary.B);
	}
	os << "b:" << endl;
	printVector(os, methodSummary.b);
	os << "x:" << endl;
	printVector(os, methodSummary.x);
	os << "error:" << endl;
	printVector(os, methodSummary.error, true);
	os << "iteration count: " << methodSummary.stepCnt << endl;
	if (methodSummary.theoreticalStepCnt != 0)
		os << "theoretical iteration count: " << methodSummary.theoreticalStepCnt << endl;
}

Vector computeError(const Matrix & A, const Vector & f, const Vector & x)
{
	Vector error(x.size(), 0), f1 = A * x;
	for (int i = 0; i < error.size(); i++)
		error[i] = abs(f1[i] - f[i]);
	return error;
}

void computeMatrixAndVectorForSimpleIteration(
	const Matrix & A, const Vector & f,
	Matrix & B, Vector & b)
{
	Matrix symmetrical = transpose(A) * A;
	double norm = computeNorm(symmetrical);
	B = loadIdentity(getN(A)) - symmetrical / norm;
	b = transpose(A) * f / norm;
}

Vector simpleIterationMethod(const Matrix & A, const Vector & f, MethodSummary & methodSummary)
{
	Matrix B;
	Vector b;	
	computeMatrixAndVectorForSimpleIteration(A, f, B, b);
	Vector x = b, oldX;

	int stepCnt = 0;
	do
	{
		oldX = x;
		x = B * oldX + b;
		stepCnt++;
	} while (!eq(x, oldX));

	methodSummary.B = B;
	methodSummary.b = b;
	methodSummary.stepCnt= stepCnt;
	methodSummary.x = x;
	methodSummary.error = computeError(A, f, methodSummary.x);
	methodSummary.theoreticalStepCnt = computeCovergenceSpeed(B, b);
	return x;
}


Vector JacobiMethod(const Matrix & A, const Vector & f, MethodSummary & methodSummary)
{
	const int n = getN(A);
	Matrix B = A;
	Vector b = f;
	for (int i = 0; i < b.size(); i++)
		b[i] /= A[i][i];
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			if (i == j)
				B[i][j] = 0;
			else
				B[i][j] = -B[i][j] / A[i][i];
	}

	int stepCnt = 0;
	Vector x = f, oldX;
	do
	{
		oldX = x;
		x = B * oldX + b;
		stepCnt++;
	} while (!eq(x, oldX));

	methodSummary.stepCnt = stepCnt;
	methodSummary.B = B;
	methodSummary.b = b;
	methodSummary.x = x;
	methodSummary.error = computeError(A, f, methodSummary.x);
	methodSummary.theoreticalStepCnt = computeCovergenceSpeed(B, b);
	return x;
}

Vector SeidelMethod(const Matrix & A, const Vector & f, MethodSummary & methodSummary)
{
	const int n = getN(A);
	Matrix B;
	Vector b;
	computeMatrixAndVectorForSimpleIteration(A, f, B, b);
	Vector x = b, oldX;
	int stepCnt = 0;
	do
	{
		oldX = x;
		for (int i = 0; i < n; i++)
		{ 
			double sum = 0;
			for (int j = 0; j < i; j++)
				sum += x[j] * B[i][j];
			for (int j = i; j < n; j++)
				sum += oldX[j] * B[i][j];
			x[i] = sum + b[i];
		}
		stepCnt++;
	} while (!eq(x, oldX));

	methodSummary.B = Matrix();
	methodSummary.b = f;
	methodSummary.error = computeError(A, f, x);
	methodSummary.stepCnt = stepCnt;
	methodSummary.theoreticalStepCnt = 0;
	methodSummary.x = x;
	return x;
}

Vector GaussSeidelMethod(const Matrix & A, const Vector & f, MethodSummary & methodSummary)
{
	const int n = getN(A);
	Vector b = f;
	for (int i = 0; i < n; i++)
		b[i] /= A[i][i];
	Matrix B = loadMatrix(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (i == j)
				B[i][j] = 0;
			else
				B[i][j] = -A[i][j] / A[i][i];
	
	Vector x = b, oldX;
	int stepCnt = 0;
	do
	{
		oldX = x;
		for (int i = 0; i < n; i++)
		{ 
			x[i] = 0;
			for (int j = 0; j < i; j++)
				x[i] += B[i][j] * x[j];
			for (int j = i; j < n; j++)
				x[i] += B[i][j] * oldX[j];
			x[i] += b[i];
		}
		stepCnt++;
	} while (!eq(x, oldX));

	methodSummary.B = B;
	methodSummary.b = b;
	methodSummary.error = computeError(A, f, x);
	methodSummary.stepCnt = stepCnt;
	methodSummary.theoreticalStepCnt = 0;
	methodSummary.x = x;
	return x;
}



int main() {
	ifstream fin("input.txt");
	ofstream fout("output.txt");
	int n;
	fin >> n;
	Matrix A = loadMatrix(n, n);
	Vector f = loadVector(n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			fin >> A[i][j];
	for (int i = 0; i < n; i++)
		fin >> f[i];

	MethodSummary methodSummary;
	simpleIterationMethod(A, f, methodSummary);
	printMethodSummary(fout, "Method of simple iteration", methodSummary);
	fout << endl;
	JacobiMethod(A, f, methodSummary);
	printMethodSummary(fout, "Jacobi method", methodSummary);
	fout << endl;
	SeidelMethod(A, f, methodSummary);
	printMethodSummary(fout, "Seidel method", methodSummary);
	fout << endl;
	GaussSeidelMethod(A, f, methodSummary);
	printMethodSummary(fout, "Gauss-Seidel method", methodSummary);
	
	return 0;
}
