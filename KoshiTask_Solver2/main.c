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
  return -(y / (x + 1) + y * y);
}

double Solution(double x)
{
  return 1 / ((x + 1) * log(x + 1));
}

vector_t FillX(double a, double h, int n)
{
  vector_t x = InitVector(n);
  int i;

  for (i = 0; i < n; ++i)
    x[i] = a + h * i;

  return x;
}

vector_t SolveKoshiTaskRungeKutta(double a, double b, double y0, vector_t x, int n, double e)
{
  vector_t y = InitVector(n);
  int i;
  double k1, k2, k3;
  double h = (b - a) / (n - 1);
  double h1 = h;
  int n1;
  int j;
  double ybuf;

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
    } while (fabs(ybuf - y[i + 1]) / 7 >= e);
    y[i + 1] = ybuf;
  }

  return y;
}

vector_t SolveKoshiExplicitAdams(double a, double b, double y0, int n, double e)
{
  vector_t y = InitVector(n);
  vector_t ystart;
  vector_t xvec;
  double h = (b - a) / (n - 1);
  double x = a;
  double k1, k2, k3;
  int i;
  
  xvec = FillX(a, h, n);
  ystart = SolveKoshiTaskRungeKutta(x, x + 2 * h, y0, xvec, 3, e);
  y[0] = ystart[0];
  y[1] = ystart[1];
  y[2] = ystart[2];
  for (i = 2; i < n - 1; ++i)
  {
    k1 = F_from_x_y(xvec[i], y[i]);
    k2 = F_from_x_y(xvec[i - 1], y[i - 1]);
    k3 = F_from_x_y(xvec[i - 2], y[i - 2]);
    y[i + 1] = y[i] + h * (23 * k1 - 16 * k2 + 5 * k3) / 12;
  }
  DeleteVector(xvec);
  DeleteVector(ystart);
  return y;
}

//vector_t SolveKoshiImplicitAdams(double a, double b, vector_t y, int n)
//{
//  double h = (b - a) / (n - 1);
//  double k1, k2, k3;
//  int i;
//  vector_t yres = InitVector(n);
//  vector_t x = FillX(a, h, n);
//
//  for (i = 0; i < n; ++i) 
//    yres[i] = y[i];
//  
//  for (i = 2; i < n; ++i)
//  {
//    k1 = F_from_x_y(x[i], yres[i]);
//    k2 = F_from_x_y(x[i - 1], yres[i - 1]);
//    k3 = F_from_x_y(x[i - 2], yres[i - 2]);
//    yres[i] = yres[i - 1] + h * (5 * k1 + 8 * k2 - k3) / 12;
//  }
//
//  DeleteVector(x);
//  return yres;
//}

vector_t SolveKoshiImplicitAdams(double a, double b, double y0, int n)
{
  double h = (b - a) / (n - 1);
  double k1, k2, k3;
  int i;
  vector_t yres = InitVector(n);
  vector_t x = FillX(a, h, n);

  vector_t ybuf = SolveKoshiTaskRungeKutta(a, a + 2 * h, y0, x, 3, 0.00001);
  yres[0] = ybuf[0];
  yres[1] = ybuf[1];
  yres[2] = ybuf[2];
  for (i = 2; i < n - 1; ++i)
  {
    k1 = F_from_x_y(x[i], yres[i]);
    k2 = F_from_x_y(x[i - 1], yres[i - 1]);
    k3 = F_from_x_y(x[i - 2], yres[i - 2]);
    yres[i + 1] = yres[i] + h * (23 * k1 - 16 * k2 + 5 * k3) / 12;

    k1 = F_from_x_y(x[i + 1], yres[i + 1]);
    k2 = F_from_x_y(x[i], yres[i]);
    k3 = F_from_x_y(x[i - 1], yres[i - 1]);
    yres[i + 1] = yres[i] + h * (5 * k1 + 8 * k2 - k3) / 12;
  }

  DeleteVector(x);
  DeleteVector(ybuf);
  return yres;
}

vector_t SolveKoshiTaskRight(vector_t x, int n)
{
  vector_t y = InitVector(n);
  int i;

  for (i = 0; i < n; ++i)
    y[i] = Solution(x[i]);

  return y;
}

double FindMax(vector_t x, int n)
{
  double max = x[0];
  int i;
  for (i = 0; i < n; ++i)
    if (max < x[i])
      max = x[i];
  return max;
}

int main(void)
{
  int n = 15;
  double a = 1;
  double b = 5;
  double h = (b - a) / (n - 1);
  double y0 = (1 / (2 * log(2))) * 1.02;
  double e = 0.00001;
  vector_t x = FillX(a, h, n);
  vector_t y1 = SolveKoshiTaskRight(x, n);
  vector_t y2 = SolveKoshiExplicitAdams(a, b, y0, n, e);
  vector_t y3 = SolveKoshiImplicitAdams(a, b, y0, n);
  int i;

  printf("Actual  \t Solved1  \tSolved2  \tError1  \tError2\n");
  for (i = 0; i < n; ++i)
  {
    printf("%g  \t%g  \t%g  \t%g  \t%g\n", y1[i], y2[i], y3[i], fabs(y1[i] - y2[i]), fabs(y1[i] - y3[i]));
  }

  {
    vector_t err1 = InitVector(n);
    vector_t err2 = InitVector(n);
    for (i = 0; i < n; ++i)
    {
      err1[i] = fabs(y1[i] - y2[i]);
      err2[i] = fabs(y1[i] - y3[i]);
    }
    printf("%g\n%g\n", FindMax(err1 + 3, n - 3), FindMax(err2 + 3, n - 3));
    DeleteVector(err1);
    DeleteVector(err2);
  }

  DeleteVector(x);
  DeleteVector(y1);
  DeleteVector(y2);
  DeleteVector(y3);

  return 0;
}