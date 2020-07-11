#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#pragma warning(disable: 4996)
#define N 9
#define A (-8.0)
#define B 0.0

int comp(const double *a, const double *b)
{
  return (*a == *b) ? 0 : ((*a < *b) ? -1 : 1);
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

int main(void)
{
  double *x;
  double *y;
  double *yder;
  int n;
  double a;
  double b;
  int i;
  double buf, buf1, max;
  double xarg;

  printf("\tEnter number of data points:\n\n");
  scanf("%i", &n);
  printf("\n\tEnter left border:\n\n");
  scanf("%lf", &a);
  printf("\n\tEnter right border:\n\n");
  scanf("%lf", &b);

  x = malloc(n * sizeof(double));
  y = malloc(n * sizeof(double));
  yder = malloc(n * sizeof(double));

  for (i = 0; i < n; i++)
  {
    x[i] = (b - a) / (n - 1) * i + a;
    y[i] = cosh(x[i]);
    yder[i] = sinh(x[i]);
//    printf("%g %g %g\n", x[i], y[i], yder[i]);
  }

  printf("\n\t\tI. Even Grid\n");
  max = 0;
  printf("\n\tCorrect polynom check (error in data points):\n\n");
  for (i = 0; i < n; i++)
  {
    buf = HermiteInterpolation(x, y, yder, n, x[i]);
    if (max < fabs(buf - y[i]))
      max = fabs(buf - y[i]);
  }
  printf("Max(eps(H(x) - f(x))) = %g\n", max);

  max = 0;
  printf("\n\tMiddle points check:\n\n");
  for (i = 0; i < n - 1; i++)
  {
    xarg = (x[i] + x[i + 1]) / 2;
    buf = HermiteInterpolation(x, y, yder, n, xarg);
    buf1 = cosh(xarg);
    if (max < fabs(buf - buf1))
      max = fabs(buf - buf1);
  }
  printf("Max(eps(H(x) - f(x))) = %g\n", max);


  printf("\n\n\t\tII. Uneven Grid\n");

  x[0] = a;
  x[1] = b;
  for (i = 2; i < n; i++)
  {
    x[i] = (double)rand() / RAND_MAX * (b - a) + a;
  }
  qsort(x, n, sizeof(double), (int (*)(const double *, const double *))comp);

  for (i = 0; i < n; i++)
  {
    y[i] = cosh(x[i]);
    yder[i] = sinh(x[i]);
//    printf("%g %g %g\n", x[i], y[i], yder[i]);
  }

  max = 0;
  printf("\n\tCorrect polynom check (error in data points):\n\n");
  for (i = 0; i < n; i++)
  {
    buf = HermiteInterpolation(x, y, yder, n, x[i]);
    if (max < fabs(buf - y[i]))
      max = fabs(buf - y[i]);
  }
  printf("Max(eps(H(x) - f(x))) = %g\n", max);

  max = 0;
  printf("\n\tMiddle points check:\n\n");
  for (i = 0; i < n - 1; i++)
  {
    xarg = (x[i] + x[i + 1]) / 2;
    buf = HermiteInterpolation(x, y, yder, n, xarg);
    buf1 = cosh(xarg);
    if (max < fabs(buf - buf1))
      max = fabs(buf - buf1);
  }
  printf("Max(eps(H(x) - f(x))) = %g\n", max);

  free(x);
  free(y);
  free(yder);

  return 0;
}