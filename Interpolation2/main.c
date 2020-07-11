#include "matrix_vector.h"
#include "slau.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#pragma warning(disable: 4996)

int SIZE;

//double Func(double x, double c)
//{
//  return cosh(x);
//}

double Func(double x, double c)
{
  return fabs(cosh(x) - cosh(c)) + cosh(c);
}

double InterFunc(vector_t coef, int m, double x)
{
  double buf = 0;
  int i;

  for (i = 0; i < m; ++i)
  {
    buf += coef[i] * exp(((((i + 2) % 2) ? 1 : -1) * (i + 1) / 2) * x);
  }

  return buf;
}

void FillMatrix(matrix_t mat, int n, int m, double a, double b, vector_t p)
{
  int i, j, k;
  double buf;
  double x;
  int i_n, j_n;
  for (i = 0; i < m; ++i)
  {
    for (j = i; j < m; ++j)
    {
      buf = 0.0;
      i_n = (((i + 2) % 2) ? 1 : -1) * (i + 1) / 2;
      j_n = (((j + 2) % 2) ? 1 : -1) * (j + 1) / 2;
      for (k = 0; k < n; ++k)
      {
        x = (b - a) / (n - 1) * k + a;
        buf += p[k] * exp(x * (i_n + j_n));
      }
      mat[i][j] = mat[j][i] = buf;
    }
  }
}

void FillVector(vector_t v, int n, int m, double a, double b, vector_t p, vector_t y)
{
  int k, i;
  double buf;
  double x;
  int i_n;

  for (i = 0; i < m; ++i)
  {
    buf = 0.0;
    i_n = (((i + 2) % 2) ? 1 : -1) * (i + 1) / 2;
    for (k = 0; k < n; ++k)
    {
      x = (b - a) / (n - 1) * k + a;
      buf += p[k] * y[k] * exp(x * i_n);
    }
    v[i] = buf;
  }
}

void FillFuncVector(vector_t y, int n, double a, double b, double c, double(*f)(double, double), double err_procent)
{
  int i;
  double x;
  double percent;

  srand(13);
  for (i = 0; i < n; ++i)
  {
    x = (b - a) / (n - 1) * i + a;
    percent = 1 + (rand() % 2 ? -1 : 1) * ((double)rand() / (RAND_MAX - 1) * err_procent);
    y[i] = f(x, c) * percent;
  }
}

int FillWeightVector(vector_t p, int m, int n, double a, double b, vector_t coef, vector_t y)
{
  int i;
  int i_max = 0;
  double max = 0;
  double x;

  for (i = 0; i < n; ++i)
  {
    x = (b - a) / (n - 1) * i + a;
    if (fabs(InterFunc(coef, m, x) - y[i]) > max)
    {
      max = fabs(InterFunc(coef, m, x) - y[i]);
      i_max = i;
    }
  }
  printf("x = %lf\n", (b - a) / (n - 1) * i_max + a);
  p[i_max] = 0.3;//(max > 1 ? 1 / (max * max) : max * max) / n;
  for (i = 0; i < i_max; ++i)
    p[i] = (1 - p[i_max]) / (n - 1);
  for (i = i_max + 1; i < n; ++i)
    p[i] = (1 - p[i_max]) / (n - 1);

  return i_max;
}

double FindErrInPoint(double(*f1)(double, double), double(*f2)(vector_t, int, double), vector_t coef, double a, double b, double c, int n, int m, int i)
{
  double x;
  x = (b - a) / (n - 1) * i + a;
  return fabs(f1(x, c) - f2(coef, m, x));
}

int main(void)
{
  vector_t f, coef1, coef2, p, y;
  matrix_t A, Q;
  int i;
  int n = 50;
  double a = -1;
  double b = 2;
  double c = 1;
  int i_change;
  int m = 5;

  //printf("  Enter number of points:\n");
  //scanf("%i", &n);
  //printf("  Enter number of basis functions:\n");
  //scanf("%i", &m);
  //printf("  Enter left border:\n");
  //scanf("%lf", &a);
  //printf("  Enter right border:\n");
  //scanf("%lf", &b);

  SetSize(m);
  f = InitVector();
  A = InitMatrix();
  SetSize(n);
  p = InitVector();
  y = InitVector();
  SetSize(m);

  FillFuncVector(y, n, a, b, c, &Func, 0.1);
  FillMatrix(A, n, m, a, b, p);
  FillVector(f, n, m, a, b, p, y);

  Q = CreateRotationMatrix(A);
  coef1 = SolveSLAU(A, Q, f);
  {
    FILE *f = fopen("out.txt", "w");
    fprintf(f, "y2 = %lf", coef1[0]);
    for (i = 1; i < m; ++i)
    {
      fprintf(f, " + (%g * exp(%i * x))", coef1[i], (((i + 2) % 2) ? 1 : -1) * (i + 1) / 2);
    }
    fclose(f);
  }

  i_change = FillWeightVector(p, m, n, a, b, coef1, y);
  FillMatrix(A, n, m, a, b, p);
  FillVector(f, n, m, a, b, p, y);
  Q = CreateRotationMatrix(A);
  coef2 = SolveSLAU(A, Q, f);

  {
    FILE *f = fopen("out.txt", "a");
    fprintf(f, "\ny3 = %lf", coef2[0]);
    for (i = 1; i < m; ++i)
    {
      fprintf(f, " + (%g * exp(%i * x))", coef2[i], (((i + 2) % 2) ? 1 : -1) * (i + 1) / 2);
    }
    fclose(f);
  }
  printf("%lf\n", FindErrInPoint(&Func, &InterFunc, coef1, a, b, c, n, m, i_change));
  printf("%lf\n", FindErrInPoint(&Func, &InterFunc, coef2, a, b, c, n, m, i_change));

  {
    FILE *f = fopen("out.txt", "a");
    double x;
    fprintf(f, "\nx1=[");
    for (i = 0; i < n; ++i)
    {
      x = (b - a) / (n - 1) * i + a;
      fprintf(f, "%lf, ", x);
    }
    fprintf(f, "];\ny5 = [");
    for (i = 0; i < n; ++i)
    {
      x = (b - a) / (n - 1) * i + a;
      fprintf(f, "%lf, ", y[i]);
    }
    fprintf(f, "];\n");
    fclose(f);
  }

  DeleteMatrix(A);
  DeleteMatrix(Q);
  DeleteVector(coef1);
  DeleteVector(coef2);
  DeleteVector(f);
  DeleteVector(p);
  DeleteVector(y);

  return 0;
}