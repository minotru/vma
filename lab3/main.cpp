//lab3, v20. Simon Karasik, course 2, group 5.
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

typedef std::vector<double> Vector;

void printVector(ostream & os, const Vector & vector,  int from, int to)
{
	for (int i = from; i < to; i++)
		os << setprecision(4) << vector[i] << '\t';
	os << endl;
}

void printVector(ostream & os, const Vector & vector)
{
	printVector(os, vector, 0, vector.size());
}

Vector solutionError(
	const Vector & a, 
	const Vector & b, 
	const Vector & c, 
	const Vector & f, 
	const Vector & y)
{
	const int n = c.size();
	Vector error(n), f1(n, 0);
	f1[0] = c[0]*y[0] - b[0]*y[1] ;
	f1[n - 1] = c[n-1]*y[n-1] - a[n-1]*y[n-2];
	for (int i = 1; i < n - 1; i++)
		f1[i] = -a[i]*y[i - 1] + c[i]*y[i] - b[i]*y[i+1];
	for (int i = 0; i < n; i++)
		error[i] = abs(f[i] - f1[i]);
	return error;
}

int main()
{
	ifstream fin("input.txt");
	ofstream fout("output.txt");


	/*
	c0 b0
	a0 c1 b1
	a1
	*/
	int n;
	fin >> n;
	Vector a(n, 0), b(n, 0), c(n, 0), f(n, 0);
	for (int i = 1; i < n; i++)
	{
		fin >> a[i];
		a[i] *= -1;
	}
	for (int i = 0; i < n - 1; i++)
	{
		fin >> b[i];
		b[i] *= -1;
	}
	for (int i = 0; i < n; i++)
		fin >> c[i];
	for (int i = 0; i < n; i++)
		fin >> f[i];
	
	Vector alpha(n + 1), beta(n + 2);

	alpha[0] = b[0] / c[0];
	beta[0] = f[0] / c[0];
	for (int i = 0; i < n; i++)
	{
		double denom = c[i] - a[i] * alpha[i];
		if (i < n - 1)
			alpha[i + 1] = b[i] / denom;
		beta[i + 1] = (f[i] + a[i] * beta[i]) / denom;
	}

	Vector y(n);
	y[n - 1] = beta[n];
	for (int i = n - 2; i >= 0; i--)
		y[i] = alpha[i + 1] * y[i + 1] + beta[i + 1];

	Vector error = solutionError(a, b, c, f, y);

	fout << "-a:" << endl;
	printVector(fout, a, 0, n - 1);
	fout << "-b:" << endl;
	printVector(fout, b, 1, n);
	fout << "c:" << endl;
	printVector(fout, c);
	fout << "alpha:" << endl;
	printVector(fout, alpha, 0, n);
	fout << "beta:" << endl;
	printVector(fout, beta, 0, n + 1);
	fout << "y:" << endl;
	printVector(fout, y);
	fout << "y error:" << endl;
	printVector(fout, error);

	return 0;
}	