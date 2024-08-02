#include "Arnoldi.hpp"

namespace lagrange {
#define MIDX(A, r, c, lda) ((A)[(r) * (lda) + (c)])

// coefficients are copied from Saad1990.pdf
#define CHEBYSHEV_DEG 14
#if CHEBYSHEV_DEG == 10

constexpr int N_ROOT = 5;

constexpr double COEF_RE[] = {0.136112052334544905E-09,
                              0.963676398167865499E+01,
                              -0.142343302081794718E+02,
                              0.513116990967461106E+01,
                              -0.545173960592769901E+00,
                              0.115698077160221179E-01};
constexpr double COEF_IM[] = {0,
                              -0.421091944767815675E+02,
                              0.176390663157379776E+02,
                              -0.243277141223876469E+01,
                              0.284234540632477550E-01,
                              0.137170141788336280E-02};
constexpr double ROOT_RE[] = {-0.402773246751880265E+01,
                              -0.328375288323169911E+01,
                              -0.171540601576881357E+01,
                              0.894404701609481378E+00,
                              0.516119127202031791E+01};
constexpr double ROOT_IM[] = {0.119385606645509767E+01,
                              0.359438677235566217E+01,
                              0.603893492548519361E+01,
                              0.858275689861307000E+01,
                              0.113751562519165076E+02};

#elif CHEBYSHEV_DEG == 14
constexpr int N_ROOT = 7;

constexpr double COEF_RE[] = {0.183216998528140087E-11,
                              0.557503973136501826E+02,
                              -0.938666838877006739E+02,
                              0.469965415550370835E+02,
                              -0.961424200626061065E+01,
                              0.752722063978321642E+00,
                              -0.188781253158648576E-01,
                              0.143086431411801849E-03};
constexpr double COEF_IM[] = {0,
                              -0.204295038779771857E+03,
                              0.912874896775456363E+02,
                              -0.116167609985818103E+02,
                              -0.264195613880262669E+01,
                              0.670367365566377770E+00,
                              -0.343696176445802414E-01,
                              0.287221133228814096E-03};
constexpr double ROOT_RE[] = {-0.562314417475317895E+01,
                              -0.508934679728216110E+01,
                              -0.399337136365302569E+01,
                              -0.226978543095856366E+01,
                              0.208756929753827868E+00,
                              0.370327340957595652E+01,
                              0.889777151877331107E+01};
constexpr double ROOT_IM[] = {0.119406921611247440E+01,
                              0.358882439228376881E+01,
                              0.600483209099604664E+01,
                              0.846173881758693369E+01,
                              0.109912615662209418E+02,
                              0.136563731924991884E+02,
                              0.166309842834712071E+02};
#endif

namespace expm {
void multiply_arnoldi_chebyshev(const std::shared_ptr<Workspace> ws,
                                size_t rate_matrix_index,
                                size_t clv_src_index,
                                size_t clv_dst_index,
                                bool transposed,
                                double t) {
  if (t == 0.0) {
#if MKL_ENABLED
    mkl_dcopy(
        ws->CLVSize(), ws->CLV(clv_src_index), 1, ws->CLV(clv_dst_index), 1);
#else
    cblas_dcopy(
        ws->CLVSize(), ws->CLV(clv_src_index), 1, ws->CLV(clv_dst_index), 1);
#endif
  }
  // allocate buffers
  constexpr int m = 20;
  const int rows = static_cast<int>(ws->matrixRows());
  const int n = rows;
  const int leading_dim = static_cast<int>(ws->leadingDimension());

  assert(rows > 0);
  assert(leading_dim > 0);

  thread_local auto A =
      std::make_unique<LagrangeMatrixBase[]>(ws->matrixSize());

  // complex buffers
  thread_local auto e1_c = std::make_unique<LagrangeComplex[]>(m);
  thread_local auto Hm_minus_theta_I_c =
      std::make_unique<LagrangeComplex[]>(m * m * N_ROOT);
  thread_local auto y_c = std::make_unique<LagrangeComplex[]>(m * N_ROOT);

  thread_local auto eHme1 =
      std::make_unique<LagrangeMatrixBase[]>(m);  // e^(Hm) * e1
  thread_local auto H = std::make_unique<LagrangeMatrixBase[]>((m + 1) * m);
  thread_local auto Q = std::make_unique<LagrangeMatrixBase[]>(n * (m + 1));
  thread_local auto e1 = std::make_unique<LagrangeMatrixBase[]>(m);

#if MKL_ENABLED
  thread_local auto ipiv_m = std::make_unique<long long int[]>(m);
#else
  thread_local auto ipiv_m = std::make_unique<int[]>(m);
#endif

  // initialization
  for (size_t i = 0; i < ws->matrixSize(); i++) {
    A.get()[i] = ws->rateMatrix(rate_matrix_index)[i] * t;
  }

  // so transpose before and after the expm operation is the same thing...
  if (transposed) {
#ifdef MKL_ENABLED
    mkl_dimatcopy('r', 't', rows, rows, 1.0, A.get(), leading_dim, leading_dim);
#else
    cblas_dimatcopy(CblasRowMajor,
                    CblasTrans,
                    rows,
                    rows,
                    1.0,
                    A.get(),
                    leading_dim,
                    leading_dim);
#endif
  }

#ifdef MKL_ENABLED
  for (size_t i = 0; i < m; ++i) {
    e1_c[i].real = 0.0;
    e1_c[i].imag = 0.0;
  }
#else
  std::fill(e1_c.get(), e1_c.get() + m, 0.0);
#endif

  lapack_complex_double_real(e1_c[0]) = 1.0;

  std::fill(e1.get(), e1.get() + m, 0.0);
  e1[0] = 1.0;

  // arnoldi iteration
  std::fill(H.get(), H.get() + (m + 1) * m, 0.0);
  std::fill(Q.get(), Q.get() + n * (m + 1), 0.0);

  const double beta = cblas_dnrm2(n, ws->CLV(clv_src_index), 1);

  cblas_dcopy(n, ws->CLV(clv_src_index), 1, Q.get(), m + 1);
  cblas_dscal(n, 1.0 / beta, Q.get(), m + 1);

  for (int j = 0; j < m; j++) {
    // candidate for the next basis: w = A * v_j+1
    cblas_dgemv(CblasRowMajor,
                CblasNoTrans,
                n,
                n,
                1.0,
                A.get(),
                n,
                Q.get() + j,
                m + 1,
                0.0,
                Q.get() + j + 1,
                m + 1);

    // remove components in directions of other bases
    for (int i = 0; i <= j; i++) {
      const double hij =
          cblas_ddot(n, Q.get() + i, m + 1, Q.get() + j + 1, m + 1);
      MIDX(H, i, j, m) = hij;
      cblas_daxpy(n, -hij, Q.get() + i, m + 1, Q.get() + j + 1, m + 1);
    }

    // compute next basis
    const double h_jp1_j = cblas_dnrm2(n, Q.get() + j + 1, m + 1);
    if (h_jp1_j > 1e-12) {
      MIDX(H, j + 1, j, m) = h_jp1_j;
      cblas_dscal(n, 1.0 / h_jp1_j, Q.get() + j + 1, m + 1);
    } else {
      break;
    }
  }

  // chebyshev approximation (see Saad1990.pdf)
  cblas_daxpby(m, COEF_RE[0], e1.get(), 1, 0.0, eHme1.get(), 1);
  const size_t m_sqrd = m * m;

  /*
   #pragma omp parallel for default(none)                            \
       shared(N_ROOT, ROOT_RE, ROOT_IM, COEF_RE, COEF_IM, m_sqrd, m, \
              Hm_minus_theta_I_c, H, y_c, e1_c, ipiv_m) num_threads(8)
  */
  for (int i = 0; i < N_ROOT; i++) {
    // note that Hm_minus_theta_I_c contains a separate buffer for each thread

    for (size_t j = 0; j < m_sqrd; j++) {
      lapack_complex_double_real(Hm_minus_theta_I_c[m_sqrd * i + j]) = -H[j];
      lapack_complex_double_imag(Hm_minus_theta_I_c[m_sqrd * i + j]) = 0.0;
    }

    LagrangeComplex theta{ROOT_RE[i], ROOT_IM[i]};
    for (size_t j = 0; j < m; j++) {
      lapack_complex_double_real(Hm_minus_theta_I_c[m_sqrd * i + j * m + j]) -=
          lapack_complex_double_real(theta);
      lapack_complex_double_imag(Hm_minus_theta_I_c[m_sqrd * i + j * m + j]) -=
          lapack_complex_double_imag(theta);
    }

    cblas_zcopy(m, e1_c.get(), 1, y_c.get() + m * i, 1);
    LAPACKE_zgesv(CblasRowMajor,
                  m,
                  1,
                  Hm_minus_theta_I_c.get() + m_sqrd * i,
                  m,
                  ipiv_m.get(),
                  y_c.get() + m * i,
                  1);

    LagrangeComplex alpha{COEF_RE[i + 1], COEF_IM[i + 1]};
    cblas_zscal(m, &alpha, y_c.get() + m * i, 1);

    // #pragma omp critical
    {
      for (size_t j = 0; j < m; j++) {
        eHme1[j] += lapack_complex_double_real(y_c[m * i + j]);
      }
    }
  }

  // update target clv directly
  cblas_dgemv(CblasRowMajor,
              CblasNoTrans,
              n,
              m,
              beta,
              Q.get(),
              m + 1,
              eHme1.get(),
              1,
              0.0,
              ws->CLV(clv_dst_index),
              1);
}

};  // namespace expm
}  // namespace lagrange
