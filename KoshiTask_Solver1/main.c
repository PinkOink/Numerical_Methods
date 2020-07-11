#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double *vector_t;

vector_t InitVector(int size)
{
  return malloc(sizeof(double) * size);
}

void DeleteVector(vector_t v)
{
  free(v);
}

double F_from_x_y(double x, double y)
{
  return -y / (x + 1) - y * y;
}

double Solution(double x)
{
  return 1 / (x + 1) / log(x + 1);
}

vector_t SolveKoshiTaskRungeKutta(double a, double b, double y0, vector_t x, int n, double e, int *nminout, int*nmaxout)
{
  vector_t y = InitVector(n);
  int i;
  double k1, k2, k3;
  double h = (b - a) / (n - 1);
  double h1 = h;
  int n1;
  int j;
  double ybuf;

  int nmax = 0;
  int nmin = INT_MAX;

  y[0] = y0;

  for (i = 0; i < n - 1; ++i)
  {
    k1 = F_from_x_y(x[i], y[i]);
    k2 = F_from_x_y(x[i] + h / 2, y[i] + h * k1 / 2);
    k3 = F_from_x_y(x[i] + h, y[i] + (2 * k2 - k1) * h);
    ybuf = y[i] + (k1 + 4 * k2 + k3) * h / 6;
    n1 = 2;
    h1 = h;
    do
    {
      h1 /= 2; 
      y[i + 1] = ybuf;
      ybuf = y[i];
      for (j = 0; j < n1; ++j)
      {
        k1 = F_from_x_y(x[i] + h1 * j, ybuf);
        k2 = F_from_x_y(x[i] + h1 * j + h1 / 2, ybuf + h1 * k1 / 2);
        k3 = F_from_x_y(x[i] + h1 * (j + 1), ybuf + (2 * k2 - k1) * h1);
        ybuf += (k1 + 4 * k2 + k3) * h1 / 6;
      }
      n1 *= 2;
    } while (fabs(ybuf - y[i + 1]) >= e);
    y[i + 1] = ybuf;

    if (n1 < nmin)
      nmin = n1;
    if (n1 > nmax)
      nmax = n1;
  }

  *nminout = nmin;
  *nmaxout = nmax;

  return y;
}

vector_t SolveKoshiTaskRight(vector_t x, int n)
{
  vector_t y = InitVector(n);
  int i;
  
  for (i = 0; i < n; ++i)
    y[i] = Solution(x[i]);

  return y;
}

vector_t FillX(double a, double h, int n)
{
  vector_t x = InitVector(n);
  int i;

  for (i = 0; i < n; ++i)
    x[i] = a + h * i;

  return x;
}

double FindMax(vector_t v, int n)
{
  double max = v[0];
  int i;
  for (i = 1; i < n; ++i)
    if (max < v[i])
      max = v[i];
  return max;
}

int main(void)
{
  int n = 21;
  double a = 1;
  double b = 5;
  double h = (b - a) / (n - 1);
  double err = 2e-4;
  double y0 = 1 / (2 * log(2)) * (1 + err);
  double e = 1e-8;
  int nmin;
  int nmax;
  vector_t x = FillX(a, h, n);
  vector_t y1 = SolveKoshiTaskRight(x, n);
  vector_t y2 = SolveKoshiTaskRungeKutta(a, b, y0, x, n, e, &nmin, &nmax);
  vector_t errv = InitVector(n);
  int i;

  for (i = 0; i < n; ++i)
  {
    printf("%g %g %g %g\n", y1[i], y2[i], fabs(y1[i] - y2[i]), 100 * fabs(y2[i] / y1[i] - 1));
    errv[i] = fabs(y1[i] - y2[i]);
  }
  printf("N max = %i\nN min = %i\n", nmax, nmin);
  printf("Error Max = %g\n", FindMax(errv, n));

  DeleteVector(x);
  DeleteVector(y1);
  DeleteVector(y2);
  DeleteVector(errv);

  return 0;
}