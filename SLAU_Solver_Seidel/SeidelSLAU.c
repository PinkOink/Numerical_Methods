#include "SeidelSLAU.h"
#include "matrix.h"
#include <stdlib.h>
#include <math.h>

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

int CheckDiagDominion(matrix_t mat)
{
  int i, j;
  double buf;

  for (i = 0; i < SIZE; i++)
  {
    buf = 0.0;
    for (j = 0; j < i; j++)
    {
      buf += fabs(mat[i][j]);
    }
    for (j = i + 1; j < SIZE; j++)
    {
      buf += fabs(mat[i][j]);
    }
    if (buf >= fabs(mat[i][i]))
    {
      return 0;
    }
  }
  
  return 1;
}

matrix_t CreateIterationMatrixD(matrix_t mat)
{
  matrix_t D = DiagMatrix(mat);
  matrix_t R = UpTriagMatrix(mat);
  matrix_t L = DownTriagMatrix(mat);

  matrix_t RL = Matrix_plus_Matrix(R, L);
  matrix_t C;
  matrix_t buf = InverseDiagMatrix(D);

  C = Matrix_x_Matrix(buf, RL);
  Matrix_x_Number(C, -1);

  DeleteMatrix(D);
  DeleteMatrix(R);
  DeleteMatrix(L);
  DeleteMatrix(RL);
  DeleteMatrix(buf);

  return C;
}

matrix_t CreateIterationMatrix(matrix_t mat, double a)
{
  matrix_t E = InitIdentityMatrix();
  matrix_t C;
  matrix_t buf = InitMatrix();

  CopyMatrix(mat, buf);
  Matrix_x_Number(buf, -a);
  C = Matrix_plus_Matrix(E, buf);

  DeleteMatrix(E);
  DeleteMatrix(buf);

  return C;
}

vector_t CreateGVector(vector_t b, double a)
{
  vector_t g = InitVector();
  
  CopyVectors(b, g);
  Vector_x_Number(g, a);

  return g;
}

vector_t CreateGVectorD(matrix_t mat, vector_t b)
{
  vector_t g;
  matrix_t D = DiagMatrix(mat);
  matrix_t D1 = InverseDiagMatrix(D);

  g = Matrix_x_Vector(D1, b);

  DeleteMatrix(D);
  DeleteMatrix(D1);

  return g;
}

void SeidelIteration(vector_t x0, matrix_t C, vector_t g, double e0, int *iterations)
{
  int i, j;
  double buf;
  vector_t xold = InitVector();
  vector_t xbuf = InitVector();

  *iterations = 0;
  do
  {
    CopyVectors(x0, xold);
    DeleteVector(xbuf);

    for (i = 0; i < SIZE; i++)
    {
      buf = 0.0;
      for (j = 0; j < SIZE; j++)
      {
        buf += x0[j] * C[i][j];
      }
      x0[i] = buf + g[i];
    }

    (*iterations)++;
    Vector_x_Number(xold, -1);
    xbuf = Vector_plus_Vector(xold, x0);
  } while (InfiniteNormVector(xbuf) >= e0);

  DeleteVector(xbuf);
  DeleteVector(xold);
}

void PrepareMatrix(matrix_t *triag)
{
  matrix_t H;
  matrix_t E = InitIdentityMatrix();
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
  H = Matrix_plus_Matrix(E, buf);             //H = I - 2 * w * w^T
  DeleteMatrix(buf);
  DeleteMatrix(E);
  buf = Matrix_x_Matrix(H, *triag);         //buf = H * D
  DeleteMatrix(*triag);
  TransposeMatrix(H);
  *triag = Matrix_x_Matrix(buf, H);          //triag = H * D * H^T
  DeleteMatrix(buf);
  DeleteMatrix(H);
}