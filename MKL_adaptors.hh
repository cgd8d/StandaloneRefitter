// Some of the MKL C interfaces are poorly designed.  Deal with that once.
#ifndef MKL_adaptors_hh
#define MKL_adaptors_hh

#include "mkl_spblas.h"
#include "mkl_types.h"
#include <limits>
#include <cassert>

// C = D*A, where D is diagonal and B is dense.
void ddiamm(size_t m, size_t n, size_t k, double* D, double* B, size_t ldb, double* C, size_t ldc)
{
  assert(m <= std::numeric_limits<MKL_INT>::max() and
         n <= std::numeric_limits<MKL_INT>::max() and
         k <= std::numeric_limits<MKL_INT>::max() and
         ldb <= std::numeric_limits<MKL_INT>::max() and
         ldc <= std::numeric_limits<MKL_INT>::max());
  char transpose = 'N';
  char* matdescra = new char[7];
  matdescra[0] = 'D';
  matdescra[1] = '*';
  matdescra[2] = 'N';
  matdescra[3] = 'F';
  matdescra[4] = '*';
  matdescra[5] = '*';
  matdescra[6] = '\0';
  MKL_INT mkl_m = m;
  MKL_INT mkl_n = n;
  MKL_INT mkl_k = k;
  MKL_INT mkl_ldb = ldb;
  MKL_INT mkl_ldc = ldc;
  double ZeroD = 0;
  double OneD = 1;
  MKL_INT ZeroI = 0;
  MKL_INT OneI = 1;

  mkl_ddiamm(&transpose, &mkl_m, &mkl_n, &mkl_k,
             &OneD, matdescra, D, &mkl_m, &ZeroI, &OneI,
             B, &mkl_ldb, &ZeroD, C, &mkl_ldc);
  delete [] matdescra;
}
#endif
