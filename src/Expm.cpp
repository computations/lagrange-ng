#include "Expm.hpp"

#include <algorithm>
#include <cmath>

#include "Common.hpp"

namespace lagrange {
#define compute_index(A, r, c, lda) ((A)[(r) * (lda) + (c)])
#define MIDX(A, r, c, lda) ((A)[(r) * (lda) + (c)])

// coefficients are copied from Saad1990.pdf
#define CHEBYSHEV_DEG 14
constexpr int N_ROOT = 7;

constexpr LagrangeComplex COEF[] = {
    {1.832169985281401e-12, 0},
    {55.75039731365018, -204.29503877977186},
    {-93.86668388770067, 91.28748967754564},
    {46.99654155503708, -11.61676099858181},
    {-9.61424200626061, -2.6419561388026267},
    {0.7527220639783216, 0.6703673655663778},
    {-0.018878125315864858, -0.03436961764458024},
    {0.00014308643141180185, 0.0002872211332288141},
};

constexpr LagrangeComplex ROOT[] = {
    {-5.623144174753179, 1.1940692161124744},
    {-5.089346797282161, 3.588824392283769},
    {-3.9933713636530257, 6.004832090996047},
    {-2.2697854309585637, 8.461738817586934},
    {0.20875692975382787, 10.991261566220942},
    {3.7032734095759565, 13.656373192499188},
    {8.897771518773311, 16.630984283471207},
};

namespace expm {

/**
 * Construct $Q$ and $H$ such that $Q  f(H)  v = f(At)v$
 *
 * Uses the normal Arnoldi process, with a modified Gram-Schmidt process to
 * normalize the basis vectors.
 *
 * For more information, please see Algorithm 6.2 in [1].
 *
 * [1]: Y. Saad, “Iterative Methods for Sparse Linear Systems”.
 *
 * @param[in] A The input matrix to factorize/reduce. A should be a square
 * matrix with dimensions `n*n` entries.
 *
 * @param[out] Q Output buffer for the $Q$ portion of the factorized matrix. $Q$
 * is a dense matrix of orthoganal basis vectors
 *
 * @param[out] H Output buffer for the $H$ portion of the factorized matrix. $H$
 * is an lower Hessenberg matrix (a mostly triangular matrix).
 *
 * @param[in] v
 * @param[in]   t
 * @param[in]   beta
 * @param[in]   n
 * @param[in]   m
 * @param[in]   leading_dim
 * @param[in]   transpose
 */
inline void arnoldi_gs(const double* A,
                       double* Q,
                       double* H,
                       const double* v,
                       double t,
                       double beta,
                       int n,
                       int m,
                       int leading_dim,
                       bool transpose) {
  constexpr double termination_epsilon = 1e-12;
  const size_t matrix_size = n * n;

  thread_local auto local_A =
      std::make_unique<LagrangeMatrixBase[]>(matrix_size);

  // initialization
  for (size_t i = 0; i < n * n; i++) { local_A[i] = A[i] * t; }

  if (transpose) {
#ifdef MKL_ENABLED
    mkl_dimatcopy('r', 't', n, n, 1.0, local_A.get(), leading_dim, leading_dim);
#else
    cblas_dimatcopy(CblasRowMajor,
                    CblasTrans,
                    n,
                    n,
                    1.0,
                    local_A.get(),
                    leading_dim,
                    leading_dim);
#endif
  }

  const int submatrix_dim = m + 1;

  // arnoldi iteration
  std::fill(Q, Q + n * submatrix_dim, 0.0);
  std::fill(H, H + submatrix_dim * m, 0.0);

  /* The first basis is compuputed outside of the loop */
  cblas_dcopy(n, v, 1, Q, submatrix_dim);
  cblas_dscal(n, 1.0 / beta, Q, submatrix_dim);

  for (int j = 0; j < m; j++) {
    double* Q_j = Q + j;
    double* w = Q + j + 1;

    /* candidate for the next basis: w = A * Q_j
     * and place w in Q_{j+1} */
    cblas_dgemv(CblasRowMajor,
                CblasNoTrans,
                n,
                n,
                1.0,
                local_A.get(),
                n,
                Q_j,
                submatrix_dim,
                0.0,
                w,
                submatrix_dim);

    // remove components in directions of other bases using Gram-Schmidt????
    for (int i = 0; i <= j; i++) {
      /* Compute the dot product between Q_i and w */
      const double hij = cblas_ddot(n, Q + i, submatrix_dim, w, submatrix_dim);

      /* Place the result in H at i, j */
      compute_index(H, i, j, m) = hij;

      /* Multiply the vector Q_i by the scalar -hij and add it to w
       * Equivalently,
       *
       *   Q_{j+1} += (Q_i, w) * Q_i
       *
       * where (.,.) is the dot product.
       *
       */
      cblas_daxpy(n, -hij, Q + i, submatrix_dim, w, submatrix_dim);
    }

    // compute next basis
    const double w_norm = cblas_dnrm2(n, w, submatrix_dim);
    /* Check if the vector is zero. If it is, the rest are also zero, and we
     * can break early. */
    if (w_norm < termination_epsilon) { break; }

    /* Basis norm goes on the diagonal */
    compute_index(H, j + 1, j, m) = w_norm;

    /* Normalize the basis vector in Q_{j+1} */
    cblas_dscal(n, 1.0 / w_norm, w, submatrix_dim);
  }
}

inline void arnoldi_householder(const double* A,
                                double* Q,
                                double* H,
                                const double* v,
                                double t,
                                double beta,
                                int n,
                                int m,
                                int leading_dim,
                                bool transpose) {
  constexpr double termination_epsilon = 1e-12;
  const size_t matrix_size = n * n;

  /* Allocate space for At */
  thread_local auto local_A =
      std::make_unique<LagrangeMatrixBase[]>(matrix_size);

  /* A := At */
  for (size_t i = 0; i < n * n; i++) { local_A[i] = A[i] * t; }

  const int submatrix_dim = m + 1;

  // arnoldi iteration
  std::fill(Q, Q + n * submatrix_dim, 0.0);
  std::fill(H, H + submatrix_dim * m, 0.0);

  /* The first basis is compuputed outside of the loop */
  cblas_dcopy(n, v, 1, Q, submatrix_dim);
  cblas_dscal(n, 1.0 / beta, Q, submatrix_dim);

  for (size_t j = 0; j < m; ++j) {
    double* Q_j = Q + j;
    double* w = Q + j + 1;

    /* candidate for the next basis: w = A * Q_j
     * and place w in Q_{j+1} */
    cblas_dgemv(CblasRowMajor,
                CblasNoTrans,
                n,
                n,
                1.0,
                local_A.get(),
                n,
                Q_j,
                submatrix_dim,
                0.0,
                w,
                submatrix_dim);

    /* Now, we need to orthogonalize Q_{j+1} */
    for (size_t i = 0; i < n; ++i) { double norm; }
  }
}

inline void chebyshev(const double* Q,
                      const double* H,
                      double* v,
                      double beta,
                      size_t n,
                      size_t m) {
  thread_local auto e1_c = std::make_unique<LagrangeComplex[]>(m);
  thread_local auto Hm_minus_theta_I_c =
      std::make_unique<LagrangeComplex[]>(m * m);
  thread_local auto y_c = std::make_unique<LagrangeComplex[]>(m);

  thread_local auto eHme1 =
      std::make_unique<LagrangeMatrixBase[]>(m);  // e^(Hm) * e1

#if MKL_ENABLED
  thread_local auto ipiv_m = std::make_unique<long long int[]>(m);
#else
  thread_local auto ipiv_m = std::make_unique<int[]>(m);
#endif

#ifdef MKL_ENABLED
  for (size_t i = 0; i < m; ++i) {
    e1_c[i].real = 0.0;
    e1_c[i].imag = 0.0;
  }
#else
  std::fill(e1_c.get(), e1_c.get() + m, 0.0);
#endif

  e1_c[0].real = 1.0;

  // chebyshev approximation (see Saad1990.pdf)
  cblas_daxpby(m, COEF[0].real, (double*)e1_c.get(), 2, 0.0, eHme1.get(), 1);

  /* This loop computes the polynomial (H - ROOT_i I) COEF */
  for (int i = 0; i < N_ROOT; i++) {
    // note that Hm_minus_theta_I_c contains a separate buffer for each thread

    /* Setup the Ith H*/
    for (size_t j = 0; j < m * m; j++) { Hm_minus_theta_I_c[j] = {-H[j], 0.0}; }

    LagrangeComplex theta = ROOT[i];
    for (size_t j = 0; j < m; j++) {
      auto val = Hm_minus_theta_I_c[j * m + j];
      Hm_minus_theta_I_c[j * m + j] = {val.real - theta.real,
                                       val.imag - theta.imag};
    }

    /* compute (H - ROOT_i * I) / e1 */
    cblas_zcopy(m, e1_c.get(), 1, y_c.get(), 1);
    LAPACKE_zgesv(CblasRowMajor,
                  m,
                  1,
                  Hm_minus_theta_I_c.get(),
                  m,
                  ipiv_m.get(),
                  y_c.get(),
                  1);

    cblas_zscal(m, COEF + i + 1, y_c.get(), 1);

    for (size_t j = 0; j < m; j++) { eHme1[j] += y_c[j].real; }
  }

  // update target clv directly
  /* Computes clv = beta * Q (e^H v)
   * Where beta = norm2(clv)
   * Using the already computed eHme1
   */
  cblas_dgemv(CblasRowMajor,
              CblasNoTrans,
              n,
              m,
              beta,
              Q,
              m + 1,
              eHme1.get(),
              1,
              0.0,
              v,
              1);
}

inline void pade(
    const double* Q, double* H, double* v, double beta, size_t n, size_t m) {
  /* Temp matrices for the computation of the exponential */
  thread_local std::unique_ptr<LagrangeMatrixBase[]> _X_1{
      new LagrangeMatrixBase[m * m]};
  thread_local std::unique_ptr<LagrangeMatrixBase[]> _X_2{
      new LagrangeMatrixBase[m * m]};
  thread_local std::unique_ptr<LagrangeMatrixBase[]> _N{
      new LagrangeMatrixBase[m * m]};
  thread_local std::unique_ptr<LagrangeMatrixBase[]> _D{
      new LagrangeMatrixBase[m * m]};
  thread_local std::unique_ptr<LagrangeMatrixBase[]> _lapack_work_buffer{
      new LagrangeMatrixBase[m]};

#ifdef MKL_ENABLED
  int ret = 0;
#endif  // MKL_ENABLED

  // We place an arbitrary limit on the size of scale exp because if it is too
  // large we run into numerical issues.
  double inf_norm = LAPACKE_dlange(CblasRowMajor, 'I', m, m, H, m);
  assert(inf_norm > 0.0);

  int scale_exp = std::min(30, std::max(0, 1 + static_cast<int>(inf_norm)));

  double Ascal = 1.0 / std::pow(2.0, scale_exp);
  assert(Ascal > 0.0);
  cblas_dscal(m * m, Ascal, H, 1);

  // q is a magic parameter that controls the number of iterations of the loop
  // higher is more accurate, with each increase of q decreasing error by 4
  // orders of magnitude. Anything above 12 is probably snake oil.
  constexpr int q = 3;
  double c = 0.5;
  double sign = -1.0;

  cblas_dcopy(m * m, H, 1, _X_1.get(), 1);
  std::fill(_X_2.get(), _X_2.get() + m * m, 0.0);
  std::fill(_N.get(), _N.get() + m * m, 0.0);
  std::fill(_D.get(), _D.get() + m * m, 0.0);

  for (size_t i = 0; i < static_cast<size_t>(m); i++) {
    _N.get()[i * m + i] = 1.0;
    _D.get()[i * m + i] = 1.0;
  }

  cblas_daxpy(m * m, c, _X_1.get(), 1, _N.get(), 1);
  cblas_daxpy(m * m, sign * c, _X_1.get(), 1, _D.get(), 1);

  // Using fortran indexing, and we started an iteration ahead to skip some
  // setup. Furthermore, we are going to unroll the loop to allow us to skip
  // some assignments.
  for (int i = 2; i <= q;) {
    c = c * (q - i + 1) / (i * (2 * q - i + 1));
    sign *= -1.0;

    cblas_dgemm(CblasRowMajor,
                CblasNoTrans,
                CblasNoTrans,
                m,
                m,
                m,
                1.0,
                H,
                m,
                _X_1.get(),
                m,
                0.0,
                _X_2.get(),
                m);
    cblas_daxpy(m * m, c, _X_2.get(), 1, _N.get(), 1);
    cblas_daxpy(m * m, sign * c, _X_2.get(), 1, _D.get(), 1);

    i += 1;

    if (i > q) { break; }

    c = c * (q - i + 1) / (i * (2 * q - i + 1));
    sign *= -1.0;

    cblas_dgemm(CblasRowMajor,
                CblasNoTrans,
                CblasNoTrans,
                m,
                m,
                m,
                1.0,
                H,
                m,
                _X_2.get(),
                m,
                0.0,
                _X_1.get(),
                m);
    cblas_daxpy(m * m, c, _X_1.get(), 1, _N.get(), 1);
    cblas_daxpy(m * m, sign * c, _X_1.get(), 1, _D.get(), 1);

    i += 1;
  }

  /* We now have _N and _D setup, so we need to compute _N * _D ^ -1*/
  {
#ifdef MKL_ENABLED
    long long int* ipiv =
        (long long int*)malloc(sizeof(long long int) * static_cast<size_t>(n));
#else
    int* ipiv = (int*)malloc(sizeof(int) * static_cast<size_t>(n));
#endif
    assert(ipiv != nullptr);

#ifdef MKL_ENABLED
    ret = LAPACKE_dgesv(CblasRowMajor, m, m, _D.get(), m, ipiv, _N.get(), m);
    assert(ret == 0);
#else
    int info = 0;
    LAPACK_dgesv(&m, &m, _D.get(), &m, ipiv, _N.get(), &m, &info);
    assert(info == 0);
#endif

    free(ipiv);
  }

  /* We now have to rescale the matrix we just computed */
  auto& r1 = _N;
  auto& r2 = _D;
  for (int i = 0; i < scale_exp; ++i) {
    cblas_dgemm(CblasRowMajor,
                CblasNoTrans,
                CblasNoTrans,
                m,
                m,
                m,
                1.0,
                r1.get(),
                m,
                r1.get(),
                m,
                0.0,
                r2.get(),
                m);
    std::swap(r1, r2);
  }

  /* r1 contains the final matrix */
  cblas_dgemv(CblasRowMajor,
              CblasNoTrans,
              n,
              m,
              beta,
              Q,
              m + 1,
              r1.get(),
              m,
              0.0,
              v,
              1);
}

void arnoldi_chebyshev(const std::shared_ptr<Workspace> ws,
                       size_t rate_matrix_index,
                       size_t clv_src_index,
                       size_t clv_dst_index,
                       bool transpose,
                       double t) {
  const auto& src_clv = ws->CLV(clv_src_index);
  auto& dst_clv = ws->CLV(clv_dst_index);

  if (t == 0.0) {
    cblas_dcopy(ws->CLVSize(), src_clv, 1, dst_clv, 1);
    return;
  }

  const auto& A = ws->rateMatrix(rate_matrix_index);

  int n = static_cast<int>(ws->matrixRows());
  constexpr int m = 20;
  int leading_dim = static_cast<int>(ws->leadingDimension());

  thread_local auto Q = std::make_unique<LagrangeMatrixBase[]>(n * (m + 1));
  thread_local auto H = std::make_unique<LagrangeMatrixBase[]>((m + 1) * m);
  const double beta = cblas_dnrm2(n, src_clv, 1);

  arnoldi_gs(
      A, Q.get(), H.get(), src_clv, t, beta, n, m, leading_dim, transpose);
  chebyshev(Q.get(), H.get(), dst_clv, beta, n, m);
}

void arnoldi_pade(const std::shared_ptr<Workspace> ws,
                  size_t rate_matrix_index,
                  size_t clv_src_index,
                  size_t clv_dst_index,
                  bool transpose,
                  double t) {
  const auto& src_clv = ws->CLV(clv_src_index);
  auto& dst_clv = ws->CLV(clv_dst_index);

  if (t == 0.0) {
    cblas_dcopy(ws->CLVSize(), src_clv, 1, dst_clv, 1);
    return;
  }

  const auto& A = ws->rateMatrix(rate_matrix_index);

  int n = static_cast<int>(ws->matrixRows());
  constexpr int m = 8;
  int leading_dim = static_cast<int>(ws->leadingDimension());

  thread_local auto Q = std::make_unique<LagrangeMatrixBase[]>(n * (m + 1));
  thread_local auto H = std::make_unique<LagrangeMatrixBase[]>((m + 1) * m);
  const double beta = cblas_dnrm2(n, src_clv, 1);

  arnoldi_gs(
      A, Q.get(), H.get(), src_clv, t, beta, n, m, leading_dim, transpose);
  pade(Q.get(), H.get(), dst_clv, beta, n, m);
}

void generate_rate_matrix(
    double* A,
    size_t lda,
    size_t restricted_state_count,
    size_t max_areas,
    size_t regions,
    std::function<double(RangeIndex, RangeIndex)> dispersion_rate,
    double extinction_rate) {
  const size_t matrix_size = lda * restricted_state_count;

  size_t source_index;
  Range source_dist;

  for (source_index = 0, source_dist = 0; source_index < restricted_state_count;
       source_dist = next_dist(source_dist, static_cast<uint32_t>(max_areas)),
      ++source_index) {
    size_t dest_index = 0;
    Range dest_dist = 0;
    for (dest_index = 0, dest_dist = 0; dest_index < restricted_state_count;
         dest_dist = next_dist(dest_dist, static_cast<uint32_t>(max_areas)),
        ++dest_index) {
      if (lagrange_popcount(source_dist ^ dest_dist) != 1) { continue; }

      if (source_dist < dest_dist && source_dist != 0) {
        double sum = 0.0;
        RangeIndex i = lagrange_fast_log2(source_dist ^ dest_dist);
        for (RangeIndex j = 0; j < regions; ++j) {
          sum += dispersion_rate(i, j) * lagrange_bextr(source_dist, j);
        }

        A[source_index * lda + dest_index] = sum;
      }

      else if (source_dist != 0) {
        A[source_index * lda + dest_index] = extinction_rate;
      }
    }
  }

  for (size_t i = 0; i < restricted_state_count; i++) {
    double sum = 0.0;
    for (size_t j = 0; j < restricted_state_count; j++) {
      sum += A[i * lda + j];
    }
    assert(sum >= 0);
    A[i * lda + i] = -sum;
  }
}

};  // namespace expm
}  // namespace lagrange
