#define _USE_MATH_DEFINES 
#include <iostream>
#include <cmath>
#include <iomanip>
#include "quanc.h"
#include "zeroin.h"
#include "fmin.h"
#include "spline.h"

double Ftemp = 0;
const short N = 21;

double zfun(double z)
{
  return std::pow(4 * z + 1, 2) - 8 * z * (std::pow(2, z) + 1) + std::pow(4, z);
}

double RungeFun(double x, double alpha)
{
  return 1.0 / (1 + alpha * x * x);
}

double integrand(double z)
{
  return 1.6 * z / (std::pow(std::sin(z), 2) + 2.56 * Ftemp * std::pow(std::cos(z), 2));
}

double FEquation(double z)
{
  Ftemp = z;
  double result = 0;
  double errest = 0;
  int nofun = 0;
  double posn = 0;
  int flag = 0;

  quanc8(integrand, 0, M_PI / 2, 1e-6, 1e-6, &result, &errest, &nofun, &posn, &flag);
  return result - 1.465862 * z;
}

double lagrange(double x, double* Xk, double* Fi)
{
	double result = 0;

	for (int i = 0; i < N; i++)
	{
	  double part = 1;
		for (int j = 0; j < N; j++)
		{
			if (i != j)
			{
				part *= (x - Xk[j]) / (Xk[i] - Xk[j]);
			}
		}
		result += part * Fi[i];
	}
	return result;
}

void calculate(double alpha, double* Xk, double* Fi)
{
	for (int i = 0; i < N; i++)
	{
		Xk[i] = -1 + i * 0.1;
		Fi[i] = RungeFun(Xk[i], alpha);
	}

	double b[N], c[N], d[N];
	int iflag = 0;
	spline(N, Xk, Fi, b, c, d, &iflag);

	std::cout << "------------------------------------------\n";
	std::cout << std::setw(10) << "X" << std::setw(10) << "f(x)" << std::setw(10) << "Lagrange" << std::setw(10) << "Spline" << std::setw(15) << "L(x) abs err" << std::setw(15) << "S(x) abs err" << "\n";
	std::cout << "------------------------------------------\n";

	for(int i = 0; i < N ; i++)
	{
		double xi = -0.95 + i * 0.1;
		double fx = RungeFun(xi, alpha);
		double Lx = lagrange(xi, Xk, Fi);
		double Sx = seval(N, xi, Xk, Fi, b, c, d, &iflag);

		std::cout << std::setw(10) << xi << std::setw(10) << fx << std::setw(10)
	            << Lx << std::setw(10) << Sx << std::setw(15) << abs(fx - Lx) << std::setw(15) << abs(fx - Sx) << "\n";
	}
	std::cout << "\n";
}

int main()
{
	double Xk[N], Fi[N];

  double Q = 80.66811 * FMin(zfun, 0.1, 0.5, 1e-6);
	std::cout << "Alpha = " << Q << "\n";
	calculate(Q, Xk, Fi);

  int flag = 0;
  double F = zeroin(0.5, 2.0, FEquation, 1e-6, &flag);
	std::cout << "Alpha = " << F << "\n";
	calculate(F, Xk, Fi);

  return 0;
}