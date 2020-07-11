#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#pragma warning(disable: 4996)

static void _swap(double *a, double *b)
{
  double buf;

  buf = *a;
  *a = *b;
  *b = buf;
}

vector_t InitVector(void)
{
  vector_t vec = calloc(SIZE, sizeof(double));
  assert(vec != NULL);
  return vec;
}

matrix_t InitMatrix(void)
{
  matrix_t mat = malloc(sizeof(double *) * SIZE);
  int i;
  assert(mat != NULL);
  for (i = 0; i < SIZE; i++)
  {
    mat[i] = malloc(sizeof(double) * SIZE);
    assert(mat[i] != NULL);
  }

  return mat;
}

matrix_t InitIdentityMatrix(void)
{
  matrix_t mat = InitMatrix();
  int i, j;

  for (i = 0; i < SIZE; i++)
  {
    for (j = 0; j < SIZE; j++)
    {
      mat[i][j] = !(i - j);
    }
  }

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

void CopyVectors(vector_t from, vector_t to)
{
  int i;

  for (i = 0; i < SIZE; i++)
  {
    to[i] = from[i];
  }
}

void CopyMatrix(matrix_t from, matrix_t to)
{
  int i, j;

  for (i = 0; i < SIZE; i++)
  {
    for (j = 0; j < SIZE; j++)
    {
      to[i][j] = from[i][j];
    }
  }
}

void PrintVector(vector_t vec)
{
  int i;

  for (i = 0; i < SIZE; i++)
  {
    printf("%lf\n", vec[i]);
  }
  printf("\n");
}

void PrintMatrix(matrix_t mat)
{
  int i, j;

  for (i = 0; i < SIZE; i++)
  {
    for (j = 0; j < SIZE; j++)
    {
      printf("%f ", mat[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

void TransposeMatrix(matrix_t mat)
{
  int i, j;

  for (i = 0; i < SIZE - 1; i++)
  {
    for (j = i + 1; j < SIZE; j++)
    {
      _swap(&mat[i][j], &mat[j][i]);
    }
  }
}

matrix_t DiagMatrix(matrix_t mat)
{
  matrix_t d = InitIdentityMatrix();
  int i;

  for (i = 0; i < SIZE; i++)
  {
    d[i][i] = mat[i][i];
  }

  return d;
}

matrix_t UpTriagMatrix(matrix_t mat)
{
  matrix_t u = InitIdentityMatrix();
  int i, j;

  for (i = 0; i < SIZE; i++)
  {
    for (j = i; j < SIZE; j++)
    {
      u[i][j] = (i == j) ? (0.0) : mat[i][j];
    }
  }

  return u;
}

matrix_t DownTriagMatrix(matrix_t mat)
{
  matrix_t d = InitIdentityMatrix();
  int i, j;

  for (i = 0; i < SIZE; i++)
  {
    for (j = 0; j <= i; j++)
    {
      d[i][j] = (i == j) ? (0.0) : mat[i][j];
    }
  }

  return d;
}

matrix_t InverseDiagMatrix(matrix_t diag)
{
  matrix_t d = InitIdentityMatrix();
  int i;

  for (i = 0; i < SIZE; i++)
  {
    d[i][i] = 1 / diag[i][i];
  }

  return d;
}

vector_t Vector_plus_Vector(vector_t A, vector_t B)
{
  vector_t C = InitVector();
  int i;

  for (i = 0; i < SIZE; i++)
  {
    C[i] = A[i] + B[i];
  }

  return C;
}

matrix_t Matrix_plus_Matrix(matrix_t A, matrix_t B)
{
  matrix_t C = InitMatrix();
  int i, j;

  for (i = 0; i < SIZE; i++)
  {
    for (j = 0; j < SIZE; j++)
    {
      C[i][j] = A[i][j] + B[i][j];
    }
  }

  return C;
}

void Vector_x_Number(vector_t vec, double n)
{
  int i;

  for (i = 0; i < SIZE; i++)
  {
    vec[i] *= n;
  }
}

void Matrix_x_Number(matrix_t mat, double n)
{
  int i, j;

  for (i = 0; i < SIZE; i++)
  {
    for (j = 0; j < SIZE; j++)
    {
      mat[i][j] *= n;
    }
  }
}

matrix_t Vector_x_VectorTRANSPOSE(vector_t vecA, vector_t vecB)
{
  matrix_t mat = InitMatrix();
  int i, j;

  for (i = 0; i < SIZE; i++)
  {
    for (j = 0; j < SIZE; j++)
    {
      mat[i][j] = vecA[i] * vecB[j];
    }
  }

  return mat;
}

vector_t Matrix_x_Vector(matrix_t mat, vector_t vec)
{
  vector_t ans = InitVector();
  int i, j;

  for (i = 0; i < SIZE; i++)
  {
    ans[i] = 0;
    for (j = 0; j < SIZE; j++)
    {
      ans[i] += mat[i][j] * vec[j];
    }
  }

  return ans;
}

matrix_t Matrix_x_Matrix(matrix_t matA, matrix_t matB)
{
  matrix_t matC = InitMatrix();
  int i, j, k;

  for (i = 0; i < SIZE; i++)
  {
    for (j = 0; j < SIZE; j++)
    {
      matC[i][j] = 0;
      for (k = 0; k < SIZE; k++)
      {
        matC[i][j] += matA[i][k] * matB[k][j];
      }
    }
  }

  return matC;
}

double EuclideanNormVector(vector_t vec)
{
  int i;
  double x = 0.0;

  for (i = 0; i < SIZE; i++)
  {
    x += vec[i] * vec[i];
  }
  return sqrt(x);
}

double InfiniteNormVector(vector_t vec)
{
  int i;
  double max = 0.0;

  for (i = 0; i < SIZE; i++)
  {
    if (fabs(vec[i]) > max)
    {
      max = fabs(vec[i]);
    }
  }

  return max;
}

double InfiniteNormMatrix(matrix_t mat)
{
  int i, j;
  double max = 0.0, buf;

  for (i = 0; i < SIZE; i++)
  {
    buf = 0.0;
    for (j = 0; j < SIZE; j++)
    {
      buf += fabs(mat[i][j]);
    }
    if (buf > max)
    {
      max = buf;
    }
  }

  return max;
}

double FirstNormMatrix(matrix_t mat)
{
  int i, j;
  double max = 0.0, buf;

  for (j = 0; j < SIZE; j++)
  {
    buf = 0.0;
    for (i = 0; i < SIZE; i++)
    {
      buf += fabs(mat[i][j]);
    }
    if (buf > max)
    {
      max = buf;
    }
  }

  return max;
}

void DeleteVector(vector_t vec)
{
  free(vec);
}

void DeleteMatrix(matrix_t mat)
{
  int i;

  for (i = 0; i < SIZE; i++)
  {
    free(mat[i]);
  }
  free(mat);
}