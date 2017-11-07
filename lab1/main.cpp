//lab1, v20. Simon Karasik, course 2, group 5.

#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace std;

typedef vector<double> Vector;
typedef vector<Vector> Matrix;

const double EPS = 10E-4;

bool eq(double x, double y) {
	return abs(x - y) < EPS;
}

Matrix loadIdentity(int n) {
	Matrix mat(n);
	for (int i = 0; i < n; i++) {
		mat[i] = Vector(n, 0);
		mat[i][i] = 1;
	}
	return mat;
}

Matrix loadMatrix(int n, int m) {
	Matrix mat(n);
	for (int i = 0; i < n; i++)
		mat[i] = Vector(m, 0);
	return mat;
}

Matrix operator*(const Matrix & a, const Matrix & b) {
	if (a[0].size() != b.size())
		throw invalid_argument("Bad size of matricies.");
	int n = a.size(), m = b[0].size(), l = a[0].size();
	Matrix p = loadMatrix(n, m);
	for (int i = 0; i < n; i++) 
		for (int j = 0; j < m; j++) {
			double sum = 0;
			for (int k = 0; k < l; k++)
				sum += a[i][k] * b[k][j];
			p[i][j] = sum;
		}	
	return p;
}

Vector operator*(const Matrix & a, const Vector & v) {
	if (a[0].size() != v.size())
		throw invalid_argument("Bad size.");
	int n = a.size(), m = v.size();
	Vector res = Vector(v.size(), 0);
	for (int i = 0; i < n; i++) 
		for (int j = 0; j < m; j++)
			res[i] += a[i][j] * v[j];
	return res;
}

void printVector(const Vector & v, ostream & os, bool useScientific = false) {
	if (useScientific) {
		os << scientific;
		for (int i = 0; i < v.size(); i++)
			os << v[i] << '\t';
	} 
	else {
		os << setprecision(4) << fixed;
		for (int j = 0; j < v.size(); j++)
			os << setw(8) << v[j];
	}
	os << endl;
}

void printMatrix(const Matrix & mat, ostream & os, bool useScientific = false) {
	for (int i = 0; i < mat.size(); i++) 
		printVector(mat[i], os, useScientific);
}

Vector Gauss(Matrix a, Vector f, double & det) {
	int n = a.size();
	if (n != f.size())
		throw invalid_argument("Size of matrix and f must be equal.");
	det = 1;

	for (int step = 0; step < n; step++) {
		int maxStep = step;
		for (int i = step + 1; i < n; i++)
			if (abs(a[i][i]) > abs(a[maxStep][maxStep]))
				maxStep = i;
		if (maxStep != step) {
			a[step].swap(a[maxStep]);
			swap(f[step], f[maxStep]);
			det *= -1;
		}
		double anchor = a[step][step];

		det *= anchor;
		f[step] /= anchor;
		for (int j = step; j < n; j++)
			a[step][j] /= anchor;
		for (int i = step + 1; i < n; i++) {
			double k = a[i][step];
			for (int j = step; j < n; j++)
				a[i][j] -= a[step][j] * k;
			f[i] -= f[step] * k;
		}

		for (int i = 0; i < n; i++) {
			Vector v = a[i];
			v.push_back(f[i]);
		}
	}

	Vector x(n, 0);
	for (int i = n - 1; i >= 0; i--) {
		x[i] = f[i];
		for (int j = n - 1; j > i; j--)
			x[i] -= a[i][j] * x[j];
	}

	return x;
}

Matrix inverse(const Matrix & a) {
	int n = a.size();
	Matrix res = loadMatrix(n, n);
	Vector f(n, 0);
	double det;
	for (int i = 0; i < n; i++) {
		f[i] = 1;
		Vector x = Gauss(a, f, det);
		for (int j = 0; j < n; j++)
			res[j][i] = x[j];
		f[i] = 0;
	}

	return res;
}

int main() {
	ifstream fin("input.txt");
	ofstream fout("output.txt");

	int n;
	fin >> n;
	Matrix a = loadMatrix(n, n);
	Vector f = Vector(n);
	for (int i = 0; i < n; i++) 
		for (int j = 0; j < n; j++)
			fin >> a[i][j];
	for (int i = 0; i < n; i++)
		fin >> f[i];
	
	double det;
	Vector x = Gauss(a, f, det);
	Matrix inversed = inverse(a);
	Vector error = a * x;
	for (int i = 0 ; i< x.size(); i++)
		error[i] -= f[i];
	Matrix inversedError = inversed * a;
	for (int i = 0; i < a.size(); i++)
		inversedError[i][i] -= 1;
	
	fout << "A:" << endl;
	printMatrix(a, fout);
	fout << "f:" << endl;
	printVector(f, fout);
	fout << "x:" << endl;
	printVector(x, fout);
	fout << "detA = " << det << endl;
	fout << "Ax - f:" << endl;
	printVector(error, fout, true);
	fout << "A^-1:" << endl;
	printMatrix(inversed, fout);
	fout << "A*A^-1 - E:" << endl;
	printMatrix(inversedError, fout, true);

	return 0;
}