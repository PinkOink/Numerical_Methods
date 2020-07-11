#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double *vector_t;
typedef vector_t *matrix_t;

vector_t InitVector(int n)
{
  return calloc(n, sizeof(double));
}

void DeleteVector(vector_t v)
{
  free(v);
}

double Solution(double x)
{
  return exp(x) / x;
}

double P(double x)
{
  return 2 / x;
}

double Q(double x)
{
  return -2;
}

double F(double x)
{
  return -exp(x) / x;
}

double U1(double x, double y, double z)
{
  return -P(x) * z - Q(x) * y;
}

double U2(double x, double y, double z)
{
  return -P(x) * z - Q(x) * y + F(x);
}

double V(double x, double y, double z)
{
  return z;
}

void PrintMatrix(matrix_t mat, int n)
{
  int i, j;
  for (i = 0; i < n; ++i)
  {
    for (j = 0; j < n; j++)
      printf("%g ", mat[i][j]);
    printf("\n");
  }
}

void PrintVector(vector_t v, int n)
{
  int i;
  for (i = 0; i < n; ++i)
    printf("%g\n", v[i]);
}

vector_t TridiagonalMatrixAlgorithm(vector_t ld, vector_t d, vector_t rd, vector_t f, int n)
{
  vector_t beta = InitVector(n);
  vector_t gamma = InitVector(n);
  vector_t x = InitVector(n);
  int i;

  beta[0] = -rd[0] / d[0];
  gamma[0] = f[0] / d[0];

  for (i = 1; i < n; ++i)
  {
    beta[i] = -rd[i] / (d[i] + ld[i] * beta[i - 1]);
    gamma[i] = (f[i] - ld[i] * gamma[i - 1]) / (d[i] + ld[i] * beta[i - 1]);
  }

  x[n - 1] = gamma[n - 1];
  for (i = n - 2; i >= 0; --i)
  {
    x[i] = beta[i] * x[i + 1] + gamma[i];
  }

  DeleteVector(beta);
  DeleteVector(gamma);

  return x;
}

vector_t SolveRight(double a, double b, int n)
{
  vector_t y = InitVector(n);
  double h = (b - a) / (n - 1);
  int i;

  for (i = 0; i < n; ++i)
    y[i] = Solution(a + h * i);

  return y;
}

vector_t SolveFiniteDifference(double a, double b, double ya, double yb, int n)
{
  double h = (b - a) / (n - 1);
  vector_t x = InitVector(n);
  vector_t f = InitVector(n);
  vector_t ld = InitVector(n);
  vector_t d = InitVector(n);
  vector_t rd = InitVector(n);
  vector_t y;
  int i;

  for (i = 0; i < n; ++i)
    x[i] = a + h * i;

  f[0] = ya;
  f[n - 1] = yb;
  for (i = 1; i < n - 1; ++i)
    f[i] = h * h * F(x[i]);

  d[0] = 1;
  d[n - 1] = 1;
  for (i = 1; i < n - 1; ++i)
  {
    ld[i] = 1 - h * P(x[i]) / 2;
    d[i] = -2 + h * h * Q(x[i]);
    rd[i] = 1 + h * P(x[i]) / 2;
  }

  y = TridiagonalMatrixAlgorithm(ld, d, rd, f, n);

  DeleteVector(x);
  DeleteVector(f);
  DeleteVector(ld);
  DeleteVector(d);
  DeleteVector(rd);

  return y;
}

vector_t SolveKoshiTaskRungeKutta(double a, double b, double y0, double z0, int n, double e, double (*U)(double, double, double), double (*V)(double, double, double))
{
  double h = (b - a) / (n - 1);
  vector_t y = InitVector(n);
  vector_t z = InitVector(n);
  double h1 = h;
  double x = a;
  int i;
  int j;
  int n1;
  double k1, k2, k3;
  double q1, q2, q3;
  double ybuf;
  double zbuf;

  y[0] = y0;
  z[0] = z0;

  for (i = 0; i < n - 1; ++i, x += h)
  {
    q1 = U(x, y[i], z[i]);
    k1 = V(x, y[i], z[i]);

    q2 = U(x + h / 2, y[i] + h * k1 / 2, z[i] + h * q1 / 2);
    k2 = V(x + h / 2, y[i] + h * k1 / 2, z[i] + h * q1 / 2);

    q3 = U(x + h, y[i] + (2 * k2 - k1) * h, z[i] + (2 * q2 - q1) * h);
    k3 = V(x + h, y[i] + (2 * k2 - k1) * h, z[i] + (2 * q2 - q1) * h);

    ybuf = y[i] + (k1 + 4 * k2 + k3) * h / 6;
    zbuf = z[i] + (q1 + 4 * q2 + q3) * h / 6;

    n1 = 2;
    h1 = h;
    do
    {
      h1 /= 2;
      y[i + 1] = ybuf;
      ybuf = y[i];
      z[i + 1] = zbuf;
      zbuf = z[i];
      for (j = 0; j < n1; ++j)
      {
        q1 = U(x + h1 * j, ybuf, zbuf);
        k1 = V(x + h1 * j, ybuf, zbuf);
        q2 = U(x + h1 * j + h1 / 2, ybuf + h1 * k1 / 2, zbuf + h1 * q1 / 2);
        k2 = V(x + h1 * j + h1 / 2, ybuf + h1 * k1 / 2, zbuf + h1 * q1 / 2);
        q3 = U(x + h1 * j + h1, ybuf + (2 * k2 - k1) * h1, zbuf + (2 * q2 - q1) * h1);
        k3 = V(x + h1 * j + h1, ybuf + (2 * k2 - k1) * h1, zbuf + (2 * q2 - q1) * h1);
        ybuf += (k1 + 4 * k2 + k3) * h1 / 6;
        zbuf += (q1 + 4 * q2 + q3) * h1 / 6;
      }
      n1 *= 2;
    } while (fabs(ybuf - y[i + 1]) / 7 >= e);
    y[i + 1] = ybuf;
    z[i + 1] = zbuf;
  }

  DeleteVector(z);

  return y;
}

vector_t Solve2KoshiTasks(double a, double b, int n, double ya, double yb, double e)
{
  double ua = 0;
  double ua2 = -1;
  double va = ya;
  double va2 = 0;
  double C;
  int i;
  vector_t y = InitVector(n);
  vector_t u = SolveKoshiTaskRungeKutta(a, b, ua, ua2, n, e, &U1, &V);
  vector_t v = SolveKoshiTaskRungeKutta(a, b, va, va2, n, e, &U2, &V);

  C = (yb - v[n - 1]) / u[n - 1];
  for (i = 0; i < n; ++i)
  {
    y[i] = C * u[i] + v[i];
  }

  DeleteVector(u);
  DeleteVector(v);

  return y;
}

int main(void)
{
  double a = 0.2;
  double b = 1;
  double ya = Solution(a);
  double yb = Solution(b);
  int n = 29;
  double x;
  double e = 1e-8;
  double h = (b - a) / (n - 1);
  vector_t y1 = SolveRight(a, b, n);
  vector_t y2 = SolveFiniteDifference(a, b, ya, yb, n);
  vector_t y3 = Solve2KoshiTasks(a, b, n, ya, yb, e);
  vector_t y4 = SolveKoshiTaskRungeKutta(a, b, ya, -20 * exp(0.2), n, e, U2, V);
  int i;

  x = a;
  for (i = 0; i < n; ++i)
  {
    printf("%g  %g  %g  %g  %g  %g\n", x, (x - a) / ((b - a) / 100), y1[i], y2[i], fabs(y1[i] - y2[i]), 100 * fabs(y2[i] / y1[i] - 1));
    x += h;
  }
  printf("\n");
  x = a;
  for (i = 0; i < n; ++i)
  {
    printf("%g  %g  %g  %g  %g  %g\n", x, (x - a) / ((b - a) / 100), y1[i], y3[i], fabs(y1[i] - y3[i]), 100 * fabs(y3[i] / y1[i] - 1));
    x += h;
  }
    printf("\n");
  x = a;
  for (i = 0; i < n; ++i)
  {
    printf("%g  %g  %g  %g  %g  %g\n", x, (x - a) / ((b - a) / 100), y1[i], y4[i], fabs(y1[i] - y4[i]), 100 * fabs(y4[i] / y1[i] - 1));
    x += h;
  }

  DeleteVector(y1);
  DeleteVector(y2);
  DeleteVector(y3);
  DeleteVector(y4);

  return 0;
}