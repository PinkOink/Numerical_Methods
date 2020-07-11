#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#pragma warning(disable: 4996)
#define N 9
#define A (-8.0)
#define B 0.0
#define pi 3.14159265358979323846

int comp(const double *a, const double *b)
{
  return (*a == *b) ? 0 : ((*a < *b) ? -1 : 1);
}

double Func(double x, double c)
{
  return cosh(x);
}

double FuncDer(double x, double c)
{
  return sinh(x);
}

double CornerFunc(double x, double c)
{
  return fabs(cosh(x) - cosh(c)) + cosh(c);
}

double CornerFuncDer(double x, double c)
{
  return (x == c) ? 0 : ((x < c) ? (-sinh(x)) : sinh(x));
}

double ChebNode(double a, double b, int i, int n)
{
  return (a + b) / 2.0 + (b - a) / 2.0 * cos(pi * (2.0 * i + 1) / 2.0 / n);
}

double AlphaJ(double x[], int n, int j)
{
  int k;
  double a = 0;

  for (k = 0; k < j; k++)
  {
    a += 1.0 / (x[j] - x[k]);
  }
  for (k = j + 1; k < n; k++)
  {
    a += 1.0 / (x[j] - x[k]);
  }

  return -2 * a;
}

double PolyFJ(double x[], int n, int j, double xarg)
{
  int i;
  double y = 1;

  for (i = 0; i < j; i++)
  {
    y *= (xarg - x[i]) / (x[j] - x[i]);
  }
  for (i = j + 1; i < n; i++)
  {
    y *= (xarg - x[i]) / (x[j] - x[i]);
  }
  return y;
}

double HermiteInterpolation(double x[], double y[], double yder[], int n, double xarg)
{
  int j;
  double herm = 0;
  double phi, psi, f;

  for (j = 0; j < n; j++)
  {
    f = PolyFJ(x, n, j, xarg);
    phi = (AlphaJ(x, n, j) * (xarg - x[j]) + 1) * f * f;
    psi = (xarg - x[j]) * f * f;

    herm += y[j] * phi + yder[j] * psi;
  }

  return herm;
}

void MakeEvenGrid(double a, double b, double c, int n, double x[], double y[], double yder[], double (*f)(double, double), double (*der)(double, double))
{
  int i;

  for (i = 0; i < n; i++)
  {
    x[i] = (b - a) / (n - 1) * i + a;
    y[i] = f(x[i], c);
    yder[i] = der(x[i], c);
  }
}

void MakeUnevenGrid(double a, double b, double c, int n, double x[], double y[], double yder[], double(*f)(double, double), double(*der)(double, double))
{
  int i;
  double alpha = 0.5;// (double)rand() / RAND_MAX;

  x[0] = a;
  x[n - 1] = b;
  for (i = 0; i < n; ++i)
  {
    x[n - 1 - i] = (1 - alpha) * ChebNode(a, b, i, n);
  }
  for (i = 0; i < n; ++i)
  {
    x[i] += alpha * ((b - a) / (n - 1) * i + a);
  }
  for (i = 0; i < n; ++i)
  {
    y[i] = f(x[i], c);
    yder[i] = der(x[i], c);
  }
}

int CoshInterpolation(int n, double a, double b)
{
  double *x;
  double *y;
  double *yder;
  double buf, buf1, max;
  int i;
  double xarg;

  printf("\t\ty = ch(x)\n");

  x = malloc(n * sizeof(double));
  y = malloc(n * sizeof(double));
  yder = malloc(n * sizeof(double));

  MakeEvenGrid(a, b, 0, n, x, y, yder, Func, FuncDer);

  printf("\n\t\tI. Even Grid\n");
  printf("\n\tCorrect polynom check (error in data points):\n\n");
  max = 0;
  for (i = 0; i < n; i++)
  {
    buf = HermiteInterpolation(x, y, yder, n, x[i]);
    if (max < fabs(buf - y[i]))
      max = fabs(buf - y[i]);
  }
  printf("Max(eps(H(x) - f(x))) = %g\n", max);

  printf("\n\tMiddle points check:\n\n");
  max = 0;
  for (i = 0; i < n - 1; i++)
  {
    xarg = (x[i] + x[i + 1]) / 2;
    buf = HermiteInterpolation(x, y, yder, n, xarg);
    buf1 = Func(xarg, 0);
 //   printf("eps = %lf\n", fabs(buf - buf1));
    if (max < fabs(buf - buf1))
      max = fabs(buf - buf1);
  }
  printf("Max(eps(H(x) - f(x))) = %g\n", max);


  printf("\n\t\tII. Uneven Grid\n");

  MakeUnevenGrid(a, b, 0, n, x, y, yder, Func, FuncDer);

  printf("\n\tCorrect polynom check (error in data points):\n\n");
  max = 0;
  for (i = 0; i < n; i++)
  {
    buf = HermiteInterpolation(x, y, yder, n, x[i]);
    if (max < fabs(buf - y[i]))
      max = fabs(buf - y[i]);
  }
  printf("Max(eps(H(x) - f(x))) = %g\n", max);

  printf("\n\tMiddle points check:\n\n");
  max = 0;
  for (i = 0; i < n - 1; i++)
  {
    xarg = (x[i] + x[i + 1]) / 2;
    buf = HermiteInterpolation(x, y, yder, n, xarg);
    buf1 = Func(xarg, 0);
//    printf("eps = %lf\n", fabs(buf - buf1));
    if (max < fabs(buf - buf1))
      max = fabs(buf - buf1);
  }
  printf("Max(eps(H(x) - f(x))) = %g\n", max);

  free(x);
  free(y);
  free(yder);

  return 0;
}

int CornerFuncInterpolation(int n, double a, double b, double c)
{
  double *x;
  double *y;
  double *yder;
  int i;
  double buf, buf1, max;
  double xarg;

  printf("\n\n\t\tCorner function\n");

  x = malloc(n * sizeof(double));
  y = malloc(n * sizeof(double));
  yder = malloc(n * sizeof(double));

  MakeEvenGrid(a, b, c, n, x, y, yder, CornerFunc, CornerFuncDer);

  printf("\n\t\tI. Even Grid\n");

  printf("\n\tCorrect polynom check (error in data points):\n\n");
  max = 0;
  for (i = 0; i < n; i++)
  {
    buf = HermiteInterpolation(x, y, yder, n, x[i]);
    if (max < fabs(buf - y[i]))
      max = fabs(buf - y[i]);
  }
  printf("Max(eps(H(x) - f(x))) = %g\n", max);

  printf("\n\tMiddle points check:\n\n");
  max = 0;
  for (i = 0; i < n - 1; i++)
  {
    xarg = (x[i] + x[i + 1]) / 2;
    buf = HermiteInterpolation(x, y, yder, n, xarg);
    buf1 = CornerFunc(xarg, c);
 //   printf("eps = %lf\n", fabs(buf - buf1));
    if (max < fabs(buf - buf1))
      max = fabs(buf - buf1);
  }
  printf("Max(eps(H(x) - f(x))) = %g\n", max);


  printf("\n\t\tII. Uneven Grid\n");

  MakeUnevenGrid(a, b, c, n, x, y, yder, CornerFunc, CornerFuncDer);

  printf("\n\tCorrect polynom check (error in data points):\n\n");
  max = 0;
  for (i = 0; i < n; i++)
  {
    buf = HermiteInterpolation(x, y, yder, n, x[i]);
    if (max < fabs(buf - y[i]))
      max = fabs(buf - y[i]);
  }
  printf("Max(eps(H(x) - f(x))) = %g\n", max);

  printf("\n\tMiddle points check:\n\n");
  max = 0;
  for (i = 0; i < n - 1; i++)
  {
    xarg = (x[i] + x[i + 1]) / 2;
    buf = HermiteInterpolation(x, y, yder, n, xarg);
    buf1 = CornerFunc(xarg, c);;
//    printf("eps = %lf\n", fabs(buf - buf1));
    if (max < fabs(buf - buf1))
      max = fabs(buf - buf1);
  }
  printf("Max(eps(H(x) - f(x))) = %g\n", max);

  free(x);
  free(y);
  free(yder);

  return 0;
}

void ResultsOut(double a, double b, double c, int n)
{
  FILE *f = fopen("out.txt", "w");
  double xarg;
  double *x;
  double *y;
  double *yder;

  x = malloc(n * sizeof(double));
  y = malloc(n * sizeof(double));
  yder = malloc(n * sizeof(double));

  fprintf(f, "x=[");
  for (xarg = a; xarg < b; xarg += 0.01)
  {
    fprintf(f, "%lf; ", xarg);
  }
  fprintf(f, "];\n\n");

  MakeEvenGrid(a, b, c, n, x, y, yder, Func, FuncDer);
  fprintf(f, "f1=[");
  for (xarg = a; xarg < b; xarg += 0.01)
  {
    fprintf(f, "%lf; ", fabs(HermiteInterpolation(x, y, yder, n, xarg) - Func(xarg, c)));
  }
  fprintf(f, "];\n\n");

  MakeUnevenGrid(a, b, c, n, x, y, yder, Func, FuncDer);
  fprintf(f, "f2=[");
  for (xarg = a; xarg < b; xarg += 0.01)
  {
    fprintf(f, "%lf; ", fabs(HermiteInterpolation(x, y, yder, n, xarg) - Func(xarg, c)));
  }
  fprintf(f, "];\n\n");


  MakeEvenGrid(a, b, c, n, x, y, yder, CornerFunc, CornerFuncDer);
  fprintf(f, "f3=[");
  for (xarg = a; xarg < b; xarg += 0.01)
  {
    fprintf(f, "%lf; ", fabs(HermiteInterpolation(x, y, yder, n, xarg) - CornerFunc(xarg, c)));
  }
  fprintf(f, "];\n\n");

  MakeUnevenGrid(a, b, c, n, x, y, yder, CornerFunc, CornerFuncDer);
  fprintf(f, "f4=[");
  for (xarg = a; xarg < b; xarg += 0.01)
  {
    fprintf(f, "%lf; ", fabs(HermiteInterpolation(x, y, yder, n, xarg) - CornerFunc(xarg, c)));
  }
  fprintf(f, "];\n\n");

  fclose(f);
  free(x);
  free(y);
  free(yder);
}

int main(void)
{
  int n = 100;
  double a = -1.0;
  double b = 2.0;
  double c = 1.0;

 /* printf("  Enter number of data points:\n\n");
  scanf("%i", &n);
  printf("\n  Enter left border:\n\n");
  scanf("%lf", &a);
  printf("\n  Enter right border:\n\n");
  scanf("%lf", &b);
  printf("\n  Enter corner point:\n\n");
  scanf("%lf", &c);*/

  CoshInterpolation(n, a, b);
  CornerFuncInterpolation(n, a, b, c);
  
  //  ResultsOut(a, b, c, n);

  return 0;
}