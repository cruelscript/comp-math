#include "spline.h"
#include <iostream>

int iflag = 0;
int last = 0;

double x[10], y[10], b[10], c[10], d[10];

double splineFunction(double arg)
{
  return seval(8, arg, x, y, b, c, d, &last) + 5 * arg - 3;
}

double(*pFunctionX)(double) = &splineFunction;

double bisection(double (*f)(double), double a, double b, double eps)
{
  while (b - a > eps)
  {
    double result = f(a) * f((b + a) / 2);

    if (result == 0) break;
    (result > 0 ? a : b) = (b + a) / 2;
  }
  return (a + b) / 2;
}
void mainb()
{
  x[0] = 0.0;
  x[1] = 0.2;
  x[2] = 0.5;
  x[3] = 0.7;
  x[4] = 1.0;
  x[5] = 1.3;
  x[6] = 1.7;
  x[7] = 2.0;

  y[0] = 1.0;
  y[1] = 1.1487;
  y[2] = 1.4142;
  y[3] = 1.6245;
  y[4] = 2.0;
  y[5] = 2.4623;
  y[6] = 3.249;
  y[7] = 4.0;

  spline(8, x, y, b, c, d, &iflag);
  std::cout << "Solution for f(x) + 5x - 3 = 0: " << bisection(pFunctionX, 0, 2, 1.0e-6);
}
