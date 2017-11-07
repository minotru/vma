//lab2, v20. Simon Karasik, course 2, group 5.
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
                sum += a[i][k] * b[k]
                       [j];
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

void printVector(const Vector & v, ostream & os, bool
                 useScientific = false) {
    if (useScientific) {
        os << scientific;
        for (int i = 0; i < v.size(); i++)
            os << v[i] << '\t';
    } else {
        os << setprecision(4) << fixed;
        for (int j = 0; j < v.size(); j++)
            os << setw(8) << v[j];
    }
    os << endl;
}
void printMatrix(const Matrix & mat, ostream & os,
                 bool useScientific = false) {
    for (int i = 0; i < mat.size(); i++)
        printVector(mat[i], os, useScientific);
}

double sqr(double x) {
    return x * x;
}

Matrix transpose(const Matrix & mat) {
    int n = mat.size(), m = mat[0].size();
    Matrix res = loadMatrix(m, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j< m; j++)
            res[j][i] = mat[i][j];
    return res;
}

int sgn(double x) {
	if (x > 0)
		return 1;
	else if (x < 0)
		return -1;
	else
		return 0;
}

Vector squareRootMethod(
	const Matrix & a, 
	const Vector & f, 
	double & det, 
	Matrix & s, 
	Vector & d, 
	Vector & y) 
{
	int n = a.size();
	s = loadIdentity(n);
    d = Vector(n);

	for (int i = 0; i < n; i++) {
		double t = 0;
		for (int k = 0; k < i; k++) 
			t += sqr(s[k][i]) * d[k];
		d[i] = sgn(a[i][i] - t);
		s[i][i] = sqrt(abs(a[i][i] - t));

		for (int j = i + 1; j < n; j++) {
			double t = 0;
			for (int k = 0; k < i; k++)
				t += s[k][i] * d[k] * s[k][j];
			s[i][j] = (a[i][j] - t) / (s[i][i] * d[i]);
		}
	}

	y = Vector(n, 0);
	for (int i = 0; i < n; i++) {
		double t = 0;
		for (int k = 0; k < i; k++)
			t += s[k][i] * y[k];
		y[i] = (f[i] - t) / s[i][i];
	}

	Vector x(n, 0);
	for (int j = n - 1; j >= 0; j--) {
		double t = 0;
		for (int k = j + 1; k < n; k++) 
			t += d[k] * s[j][k] * x[k];
		x[j] = (y[j] - t) / (d[j] * s[j][j]);
	}
	
	det = 1;
	for (int i = 0; i < n; i++)
		det *= d[i] * sqr(s[i][i]);

	return x;
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
	f = transpose(a) * f;
    a = transpose(a) * a;
	double det;
	Matrix s;
	Vector d, y;
	const Vector x = squareRootMethod(a, f, det, s, d, y);	
    det = sqrt(det);
    Vector error = a * x;
    for (int i = 0 ; i< x.size(); i++)
        error[i] -= f[i];
    
	fout << "A^T * A:" << endl;
	printMatrix(a, fout);
    fout << "A^T * f:" << endl;
    printVector(f, fout);
	fout << "S:" << endl;
	printMatrix(s, fout);
	fout << "diag(D):" << endl;
	printVector(d, fout);
	fout << "y:" << endl;
	printVector(y, fout);
    fout << "x:" << endl;
    printVector(x, fout);
    fout << "|detA| = " << det << endl;
    fout << "error of x:" << endl;
    printVector(error, fout, true);
    
	return 0;
}
