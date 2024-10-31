#include "Expm.hpp"


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
 * Construct Q and H such that Q * f(H) * v = f(At) * v
 */
void arnoldi(const double* A,
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
    // candidate for the next basis: w = A * q_j+1
    cblas_dgemv(CblasRowMajor,
                CblasNoTrans,
                n,
                n,
                1.0,
                local_A.get(),
                n,
                Q + j,
                submatrix_dim,
                0.0,
                Q + j + 1,
                submatrix_dim);

    // remove components in directions of other bases
    for (int i = 0; i <= j; i++) {
      /* Compute the dot product between Q_i and Q_j */
      const double hij =
          cblas_ddot(n, Q + i, submatrix_dim, Q + j + 1, submatrix_dim);

      /* Place the result in H at i, j */
      compute_index(H, i, j, m) = hij;

      /* Multiply the vector Q_i by the scalar -hij and add it to Q_{j+1} */
      cblas_daxpy(n, -hij, Q + i, submatrix_dim, Q + j + 1, submatrix_dim);
    }

    // compute next basis
    const double Q_j1_norm = cblas_dnrm2(n, Q + j + 1, submatrix_dim);
    /* Check if the vector is zero. If it is, the rest are also zero, and we
     * can break early. */
    if (Q_j1_norm < termination_epsilon) { break; }

    /* Basis norm goes on the diagonal */
    compute_index(H, j + 1, j, m) = Q_j1_norm;

    /* Normalize the basis vector in Q_{j+1} */
    cblas_dscal(n, 1.0 / Q_j1_norm, Q + j + 1, submatrix_dim);
  }
}

void chebyshev(const double* Q,
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

void arnoldi_chebyshev(const std::shared_ptr<Workspace> ws,
                     size_t rate_matrix_index,
                     size_t clv_src_index,
                     size_t clv_dst_index,
                     bool transpose,
                     double t) {
  const auto& A = ws->rateMatrix(rate_matrix_index);
  const auto& src_clv = ws->CLV(clv_src_index);
  auto& dst_clv = ws->CLV(clv_dst_index);

  int n = static_cast<int>(ws->matrixRows());
  constexpr int m = 20;
  int leading_dim = static_cast<int>(ws->leadingDimension());

  thread_local auto Q = std::make_unique<LagrangeMatrixBase[]>(n * (m + 1));
  thread_local auto H = std::make_unique<LagrangeMatrixBase[]>((m + 1) * m);
  const double beta = cblas_dnrm2(n, src_clv, 1);

  arnoldi(A, Q.get(), H.get(), src_clv, t, beta, n, m, leading_dim, transpose);
  chebyshev(Q.get(), H.get(), dst_clv, beta, n, m);
}

};  // namespace expm
}  // namespace lagrange