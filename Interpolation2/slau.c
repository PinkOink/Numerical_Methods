#include "slau.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#pragma warning(disable:4996)

extern int SIZE;

matrix_t InitTriagMatrix(void)
{
  matrix_t mat = InitIdentityMatrix();
  int i, j;

  for (i = 0; i < SIZE; i++)
  {
    for (j = i; j < SIZE; j++)
    {
      mat[i][j] = (double)rand() / RAND_MAX * 10000 + 1;
    }
  }
 // PrintMatrix(mat);
  return mat;
}

matrix_t GetMatrix(FILE *f)
{
  matrix_t mat = InitIdentityMatrix();
  int i, j;

  for (i = 0; i < SIZE; i++)
  {
    for (j = 0; j < SIZE; j++)
    {
      fscanf(f, "%lf", &mat[i][j]);
    }
  }

  return mat;
}

vector_t InitRandVector(void)
{
  vector_t vec = InitVector();
  int i;

  for (i = 0; i < SIZE; i++)
  {
    vec[i] = (double)rand() / RAND_MAX * 100 + 1;
  }

  return vec;
}

void PrepareMatrix(matrix_t *diag)
{
  matrix_t H = InitIdentityMatrix();
  matrix_t buf;
  vector_t w = InitVector();
  double n = 1.0;
  int i;

  for (i = 0; i < SIZE - 1; i++)
  {
    w[i] = ((double)rand() / RAND_MAX) * n;
    n -= w[i] * w[i];
  }
  w[i] = sqrt(n);

//  printf("\n%lf\n", EuclideanNormVector(w));

  buf = Vector_x_VectorTRANSPOSE(w, w);   //buf = w * w^T  
  DeleteVector(w);
  Matrix_x_Number(buf, -2.0);             //buf = -2 * w * w^T
  Matrix_plus_Matrix(H, buf);             //H = I - 2 * w * w^T
  DeleteMatrix(buf);                      
  buf = Matrix_x_Matrix(H, *diag);         //buf = H * D
  DeleteMatrix(*diag);
  TransposeMatrix(H);
  *diag = Matrix_x_Matrix(buf, H);          //diag = H * D * H^T
  DeleteMatrix(buf);
  DeleteMatrix(H);
}

vector_t DeltaVector(vector_t vec)
{
  vector_t delta = InitVector();
  int i;

  for (i = 0; i < SIZE; i++)
  {
    delta[i] = ((double)rand() / RAND_MAX) * vec[i] / 100.0;
  }

  return delta;
}

matrix_t DeltaMatrix(matrix_t mat)
{
  matrix_t delta = InitMatrix();
  int i, j;

  for (i = 0; i < SIZE; i++)
  {
    for (j = 0; j < SIZE; j++)
    {
      delta[i][j] = ((double)rand() / RAND_MAX) * mat[i][j] / 100.0;
    }
  }

  return delta;
}

matrix_t InitRotationMatrix(int i, int j, double tg)
{
  matrix_t mat = InitIdentityMatrix();
  double cos = 1.0 / sqrt(1 + tg * tg);
  double sin = tg * cos;

  mat[i][i] = cos;
  mat[i][j] = sin;
  mat[j][i] = -sin;
  mat[j][j] = cos;

  return mat;
}

matrix_t CreateRotationMatrix(matrix_t mat)
{
  int i, j;
  matrix_t Q = InitIdentityMatrix();
  matrix_t Q1;
  matrix_t A = InitIdentityMatrix();
  matrix_t buf1 = Matrix_x_Matrix(A, mat);
  matrix_t buf2;
  DeleteMatrix(A);
  A = buf1;

  for (i = 0; i < SIZE - 1; i++)
  {
    for (j = i + 1; j < SIZE; j++)
    {
      Q1 = InitRotationMatrix(i, j, A[j][i] / A[i][i]);
      buf1 = Matrix_x_Matrix(Q1, Q);
      buf2 = Matrix_x_Matrix(Q1, A);
      DeleteMatrix(Q);
      DeleteMatrix(Q1);
      DeleteMatrix(A);
      Q = buf1;
      A = buf2;
    }
  }
  DeleteMatrix(A);

  return Q;
}

vector_t SolveSLAU(matrix_t mat, matrix_t rotat, vector_t b)
{
  vector_t x = InitVector();
  matrix_t A1 = Matrix_x_Matrix(rotat, mat);
  vector_t B1 = Matrix_x_Vector(rotat, b);
  int i, j;
  
  for (i = SIZE - 1; i >= 0; i--)
  {
    x[i] = B1[i];

    for (j = SIZE - 1; j > i; j--)
    {
      x[i] -= A1[i][j] * x[j];
    }

    x[i] /= A1[i][i];
  }

  DeleteMatrix(A1);
  DeleteVector(B1);

  return x;
}