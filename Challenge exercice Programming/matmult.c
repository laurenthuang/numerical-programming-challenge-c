#include "matmult.h"
#define A_lda(i, j) A[(i) * K + (j)]
#define B_lda(i, j) B[(i) * N + (j)]
#define C_lda(i, j) C[(i) * N + (j)]

#define min(x, y) ((x) < (y) ? (x) : (y))

void matmult(const int64_t M, const int64_t N, const int64_t K,
             const double *const A, const double *const B, double *const C)
{
  // from test 2. i discovered BLOCK_SIZE = 16 or 32 are faster, especially 16.
  const int64_t BLOCK_SIZE = 16;
  double B_block[BLOCK_SIZE][BLOCK_SIZE];
  for (int64_t i = 0; i < M; i++)
  {
    for (int64_t j = 0; j < N; j++)
    {
      C_lda(i, j) = 0.0; // initialize
    }
  }
  for (int ii = 0; ii < M; ii += BLOCK_SIZE)
  {
    int ib = min(M - ii, BLOCK_SIZE);
    for (int jj = 0; jj < N; jj += BLOCK_SIZE)
    {
      int jb = min(N - jj, BLOCK_SIZE);
      for (int pp = 0; pp < K; pp += BLOCK_SIZE)
      {
        int pb = min(K - pp, BLOCK_SIZE);
        for (int j = 0; j < jb; j++)
        {
          for (int p = 0; p < pb; p++)
          {
            B_block[j][p] = B_lda(pp + p, jj + j); // B transpose
          }
        }
        for (int i = 0; i < ib; i++)
        {
          for (int j = 0; j < jb; j++)
          {
            for (int p = 0; p <= pb - 4; p += 4) // try divid the loop since 4 is the LCD
            {
              C_lda(ii + i, jj + j) += A_lda(ii + i, pp + p) * B_block[j][p] + A_lda(ii + i, pp + p + 1) * B_block[j][p + 1] + A_lda(ii + i, pp + p + 2) * B_block[j][p + 2] + A_lda(ii + i, pp + p + 3) * B_block[j][p + 3];
            }
          }
        }
      }
    }
  }
}