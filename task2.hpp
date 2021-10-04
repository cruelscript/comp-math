#include <iostream>
#include <iomanip>
#include "decomp.h"

const size_t N = 5;

int ipvt[N];
double A[N][N];
double AAR[N][N];
double R[N][N];
double a[N - 1] = {4, 3, 2, 1.5};
double cond = 0;
int flag = 0;
double norm = 0;

void createMatrix(double(&A)[N][N], double(&a)[N - 1])
{
  for (size_t i = 0; i < N; ++i)
  {
    for (size_t j = 0; j < N; ++j)
    {
      if (i >= j)
      {
        A[i][j] = 1;
      }
      else
      {
        A[i][j] = a[i];
      }
    }
  }
}

void createEMatrix(double(&A)[N][N])
{
  for (size_t i = 0; i < N; ++i)
  {
    for (size_t j = 0; j < N; ++j)
    {
      if (i == j)
      {
        A[i][j] = 1;
      }
      else
      {
        A[i][j] = 0;
      }
    }
  }
}

void reversedMatrix(double(&A)[N][N], double(&AR)[N][N])
{
  double E[N][N];
  double b[N];

  createEMatrix(E);

  for (size_t i = 0; i < N; ++i)
  {
    for (size_t j = 0; j < N; ++j)
    {
      b[j] = E[i][j];
    }
    solve(5, 5, *A, b, ipvt);
    for (size_t j = 0; j < N; ++j)
    {
      AR[j][i] = b[j];
    }
  }
}

void multiplyMatrixes(double(&A)[N][N], double(&B)[N][N], double(&result)[N][N])
{
  double sum = 0;
  for (size_t i = 0; i < N; ++i)
  {
    for (size_t j = 0; j < N; ++j)
    {
      sum = 0;
      for (size_t k = 0; k < N; ++k)
      {
        sum += A[i][k] * B[k][j];
        result[i][j] = sum;
      }
    }
  }
}

void findRMatrix(double(&AAR)[N][N], double(&R)[N][N], double& norm)
{
  double E[N][N];
  createEMatrix(E);

  for (size_t i = 0; i < N; ++i)
  {
    for (size_t j = 0; j < N; ++j)
    {
      R[i][j] = AAR[i][j] - E[i][j];
    }
  }

  double sum = 0;
  norm = 0;
  for (size_t i = 0; i < N; ++i)
  {
    sum = 0;
    for (size_t j = 0; j < N; ++j)
    {
      sum += abs(R[i][j]);
    }
    norm = norm < sum ? sum : norm;
  }
}

void printMatrix(double(&A)[N][N])
{
  for (size_t i = 0; i < N; ++i)
  {
    for (size_t j = 0; j < N; ++j)
    {
      std::cout << std::setprecision(5) << std::setiosflags(std::ios::scientific) << A[i][j] << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}

void computeTask2(double a4)
{
  a[3] = a4;
  std::cout << "Matrix A(" << a4 << "):\n";
  createMatrix(A, a);
  printMatrix(A);

  double Aold[N][N];
  for (size_t i = 0; i < N; ++i)
  {
    for (size_t j = 0; j < N; ++j)
    {
      Aold[i][j] = A[i][j];
    }
  }

  decomp(N, N, *A, &cond, ipvt, &flag);
  double AR[N][N];
  reversedMatrix(A, AR);
  std::cout << "Matrix A^(-1):\n";
  printMatrix(AR);

  multiplyMatrixes(Aold, AR, AAR);
  findRMatrix(AAR, R, norm);
  std::cout << "Matrix R:\n";
  printMatrix(R);
  std::cout << "Cond = " << cond << " Norm = " << norm << "\n\n";
  std::cout << "______________________________________________\n";
}

int nemain()
{
  computeTask2(1.5);
  computeTask2(1.01);
  computeTask2(1.001);
  computeTask2(1.0001);

  return 0;
}