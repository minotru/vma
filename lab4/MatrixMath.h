#pragma once

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace std;

typedef vector<double> Vector;
typedef vector<Vector> Matrix;

int getN(const Matrix &);
int getM(const Matrix &);

Vector loadVector(int n)
{
	return Vector(n, 0);
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

Matrix operator*(Matrix m, double v)
{
	for (int i = 0; i < getN(m); i++)
		for (int j = 0; j < getM(m); j++)
			m[i][j] *= v;
	return m;
}

Vector operator*(Vector vec, double v)
{
	for (int i = 0; i < vec.size(); i++)
		vec[i] *= v;
	return vec;
}

Vector operator+(const Vector & v1, const Vector & v2)
{
	Vector res = v1;
	for (int i = 0; i < v1.size(); i++)
		res[i] += v2[i];
	return  res;
}

Vector operator/(const Vector &  vec, double v)
{
	return vec * (1 / v);
}

Matrix operator*(double v, const Matrix & m)
{
	return m * v;
}

Matrix operator/(const Matrix & m, double v)
{
	return m * (1 / v);
}

Matrix operator-(Matrix m)
{
	return m * (-1);
}

Matrix operator+(const Matrix & a, const Matrix & b)
{
	if (getN(a) != getN(b) || getM(a) != getM(b))
		throw invalid_argument("size(a) != size(b), can not sum");
	Matrix res(a);
	for (int i = 0; i < getN(b); i++)
		for (int j = 0; j < getM(b); j++)
			res[i][j] += b[i][j];
	return res;
}

Matrix operator-(const Matrix & a, const Matrix & b)
{
	return a + (-b);
}


void printVector(ostream & os, const Vector & v, bool useScientific = false)
{
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
void printMatrix(ostream & os, const Matrix & mat, bool useScientific = false) 
{
	for (int i = 0; i < mat.size(); i++)
		printVector(os,  mat[i], useScientific);
}

Matrix transpose(const Matrix & mat) {
	int n = mat.size(), m = mat[0].size();
	Matrix res = loadMatrix(m, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j< m; j++)
			res[j][i] = mat[i][j];
	return res;
}

int getN(const Matrix & mat)
{
	return mat.size();
}

int getM(const Matrix & mat)
{
	return mat[0].size();
}