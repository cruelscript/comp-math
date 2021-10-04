#include <iostream>
#include "rkf45.h"
#include <cmath>
#include <iomanip>

int myDiffEq(int NEQN, double t, double* x, double* dx)
{
  dx[0] = -73 * x[0] - 210 * x[1] + log(1 + t * t);
  dx[1] = x[0] + exp(-t) + t * t + 1;
  return NEQN;
}

void EulerCauchy(double t, double* zn, double* zn1, double h)
{
  double k1[2];
  double k2[2];
  double tempk[2];

  myDiffEq(2, t, zn, k1);
  k1[0] *= h;
  k1[1] *= h;

  tempk[0] = zn[0] + k1[0];
  tempk[1] = zn[1] + k1[1];
  myDiffEq(2, t + h, tempk, k2);
  k2[0] *= h;
  k2[1] *= h;

  zn1[0] = zn[0] + (k1[0] + k2[0]) / 2;
  zn1[1] = zn[1] + (k1[1] + k2[1]) / 2;
}
int nemainb()
{
  double y0[] = {-3, 1};
  double yp[2];

  double h, relerr, abserr, x1, x2;
  int n, flag, nfe, maxfe, fail, step;
  n = 2;
  flag = 1;
  maxfe = 5000;
  relerr = 1.0e-4;
  abserr = 1.0e-4;

  rkfinit(n, &fail);

  if (fail == 0)
  {
    std::cout << "   t      x1           x2\n";
    std::cout << "---------------------------------\n";

    for (step = 1; step <= 1 / 0.05; ++step)
    {
      x2 = 0.05 * step;
      x1 = x2 - 0.05;
      rkf45(myDiffEq, n, y0, yp, &x1, x2, &relerr, abserr, &h, &nfe, maxfe, &flag);
      std::cout << std::setw(5) << std::setprecision(8) << x2 << std::setw(14) << y0[0] << std::setw(14) << y0[1] << "\n";
    }
  }
  rkfend();
  std::cout << "\n" << nfe << "\n" << h;

  std::cout << "\n\nEuler-Cauchy, h = 0.05\n";

  double z1[] = {-3, 1};
  double z2[2];

  std::cout << "   t      x1           x2\n";
  std::cout << "---------------------------------\n";

  double st = 0;
  while (st <= 1.01)
  {
    EulerCauchy(st, z1, z2, 0.05);
    std::cout << std::setw(5) << std::setprecision(8) << st << std::setw(14) << z1[0] << std::setw(14) << z1[1] << "\n";
    st += 0.05;
    z1[0] = z2[0];
    z1[1] = z2[1];
  }

  std::cout << "\n\nEuler-Cauchy, h = 0.027\n";

  z1[0] = -3;
  z1[1] = 1;


  std::cout << "   t      x1           x2\n";
  std::cout << "---------------------------------\n";

  st = 0;
  while (st <= 1.01)
  {
    EulerCauchy(st, z1, z2, 0.0275);
    std::cout << std::setw(5) << std::setprecision(8) << st << std::setw(14) << z1[0] << std::setw(14) << z1[1] << "\n";
    st += 0.05;
    z1[0] = z2[0];
    z1[1] = z2[1];
  }

  return 0;
}
