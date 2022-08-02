//
// Created by sunflower on 24.07.22.
//

#ifndef LAGRANGE_CPP_ARNOLDI_H
#define LAGRANGE_CPP_ARNOLDI_H

#include <algorithm>
#include <memory>

#include "Common.h"
#include "Workspace.h"

/*
 * expm_multiply_arnoldi_chebyshev:
 *   combination of arnoldi and chebyshev approximation
 *
 * expm_multiply_arnoldi_pade:
 *   combination of arnoldi and pade (scaling and squaring)
 *
 * expm_multiply_pade:
 *   pade (scaling and squaring)
 *
 * Note that all three functions can be directly substituted in Operation.cpp,
 * line 544.
 *
 * Referenced papers: Saad1990.pdf and Saad1992.pdf
 */

namespace arnoldi {

// coefficients are copied from Saad1990.pdf
// clang-format off
#define CHEBYSHEV_DEG 14
#if CHEBYSHEV_DEG == 10
constexpr double COEF_RE[] = {0.136112052334544905E-09, 0.963676398167865499E+01, -0.142343302081794718E+02, 0.513116990967461106E+01, -0.545173960592769901E+00, 0.115698077160221179E-01};
constexpr double COEF_IM[] = {0, -0.421091944767815675E+02, 0.176390663157379776E+02, -0.243277141223876469E+01, 0.284234540632477550E-01, 0.137170141788336280E-02};
constexpr double ROOT_RE[] = {-0.402773246751880265E+01, -0.328375288323169911E+01, -0.171540601576881357E+01, 0.894404701609481378E+00, 0.516119127202031791E+01};
constexpr double ROOT_IM[] = {0.119385606645509767E+01, 0.359438677235566217E+01, 0.603893492548519361E+01, 0.858275689861307000E+01, 0.113751562519165076E+02};
constexpr int N_ROOT = 5;
constexpr int N_COEF = 6;
#elif CHEBYSHEV_DEG == 14
constexpr double COEF_RE[] = {0.183216998528140087E-11, 0.557503973136501826E+02, -0.938666838877006739E+02, 0.469965415550370835E+02, -0.961424200626061065E+01, 0.752722063978321642E+00, -0.188781253158648576E-01, 0.143086431411801849E-03};
constexpr double COEF_IM[] = {0, -0.204295038779771857E+03, 0.912874896775456363E+02, -0.116167609985818103E+02, -0.264195613880262669E+01, 0.670367365566377770E+00, -0.343696176445802414E-01, 0.287221133228814096E-03};
constexpr double ROOT_RE[] = {-0.562314417475317895E+01, -0.508934679728216110E+01, -0.399337136365302569E+01, -0.226978543095856366E+01, 0.208756929753827868E+00, 0.370327340957595652E+01, 0.889777151877331107E+01};
constexpr double ROOT_IM[] = {0.119406921611247440E+01, 0.358882439228376881E+01, 0.600483209099604664E+01, 0.846173881758693369E+01, 0.109912615662209418E+02, 0.136563731924991884E+02, 0.166309842834712071E+02};
constexpr int N_ROOT = 7;
constexpr int N_COEF = 8;
#endif
// clang-format on

#define MIDX(A, r, c, lda) ((A)[(r) * (lda) + (c)])

// see Saad1990.pdf and Saad1992.pdf
void expm_multiply_arnoldi_chebyshev(const std::shared_ptr<Workspace> ws,
                                     size_t rate_matrix, size_t clv_src,
                                     size_t clv_dst, bool transposed,
                                     size_t prob_matrix, double t) {
  //
  // allocate buffers
  //

  const int m = 20;
  const int rows = static_cast<int>(ws->matrix_rows());
  const int n = rows;
  int leading_dim = static_cast<int>(ws->leading_dimension());

  assert(rows > 0);
  assert(leading_dim > 0);

  auto A = std::make_unique<lagrange_matrix_base_t[]>(ws->matrix_size());

  // complex buffers
  auto e1_c = std::make_unique<lagrange_complex_t[]>(m);
  auto Hm_minus_theta_I_c =
      std::make_unique<lagrange_complex_t[]>(m * m * N_ROOT);
  auto y_c = std::make_unique<lagrange_complex_t[]>(m * N_ROOT);

  auto eHme1 = std::make_unique<lagrange_matrix_base_t[]>(m);  // e^(Hm) * e1
  auto H = std::make_unique<lagrange_matrix_base_t[]>((m + 1) * m);
  auto Q = std::make_unique<lagrange_matrix_base_t[]>(n * (m + 1));
  auto e1 = std::make_unique<lagrange_matrix_base_t[]>(m);

  auto ipiv_m = std::make_unique<int[]>(n);

  //
  // initialization
  //

  for (size_t i = 0; i < ws->matrix_size(); i++) {
    A.get()[i] = ws->rate_matrix(rate_matrix)[i] * t;
  }

  // so transpose before and after the expm operation is the same thing...
  if (transposed) {
    cblas_dimatcopy(CblasRowMajor, CblasTrans, rows, rows, 1.0, A.get(),
                    leading_dim, leading_dim);
  }

  std::fill(e1_c.get(), e1_c.get() + m, 0.0);
  e1_c[0] = 1.0;

  std::fill(e1.get(), e1.get() + m, 0.0);
  e1[0] = 1.0;

  //
  // arnoldi iteration
  //

  std::fill(H.get(), H.get() + (m + 1) * m, 0.0);
  std::fill(Q.get(), Q.get() + n * (m + 1), 0.0);

  const double beta = cblas_dnrm2(n, ws->clv(clv_src), 1);

  cblas_dcopy(n, ws->clv(clv_src), 1, Q.get(), m + 1);
  cblas_dscal(n, 1.0 / beta, Q.get(), m + 1);

  for (int j = 0; j < m; j++) {
    // candidate for the next basis: w = A * v_j+1
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, A.get(), n, Q.get() + j,
                m + 1, 0.0, Q.get() + j + 1, m + 1);

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

  // Sorry it seems that I forgot to test with OpenMP...
  // #pragma omp parallel for default(none)                            \
//     shared(N_ROOT, ROOT_RE, ROOT_IM, COEF_RE, COEF_IM, m_sqrd, m, \
//            Hm_minus_theta_I_c, H, y_c, e1_c, ipiv_m) num_threads(8)
  for (int i = 0; i < N_ROOT; i++) {

    // note that Hm_minus_theta_I_c contains a separate buffer for each thread
    for (size_t j = 0; j < m_sqrd; j++)
      Hm_minus_theta_I_c[m_sqrd * i + j] = -H[j];

    lagrange_complex_t theta(ROOT_RE[i], ROOT_IM[i]);
    for (size_t j = 0; j < m; j++) {
      Hm_minus_theta_I_c[m_sqrd * i + j * m + j] -= theta;
    }

    cblas_zcopy(m, e1_c.get(), 1, y_c.get() + m * i, 1);
    LAPACKE_zgesv(CblasRowMajor, m, 1, Hm_minus_theta_I_c.get() + m_sqrd * i, m,
                  ipiv_m.get(), y_c.get() + m * i, 1);

    lagrange_complex_t alpha(COEF_RE[i + 1], COEF_IM[i + 1]);
    cblas_zscal(m, &alpha, y_c.get() + m * i, 1);

    // #pragma omp critical
    {
      for (size_t j = 0; j < m; j++)
        eHme1[j] += lapack_complex_double_real(y_c[m * i + j]);
    }
  }

  // update target clv directly
  cblas_dgemv(CblasRowMajor, CblasNoTrans, n, m, beta, Q.get(), m + 1,
              eHme1.get(), 1, 0.0, ws->clv(clv_dst), 1);
}

// old good Pade method (this is used in expm_multiply_arnoldi_pade
// and expm_multiply_pade
auto expm_pade(const std::unique_ptr<lagrange_matrix_base_t[]> &A, size_t n,
               size_t leading_dim, bool transposed, double t)
    -> std::unique_ptr<lagrange_matrix_base_t[]> {
  int ret = 0;
  int rows = n;

  assert(rows > 0);

  auto X_1 = std::make_unique<lagrange_matrix_base_t[]>(n * n);
  auto X_2 = std::make_unique<lagrange_matrix_base_t[]>(n * n);
  auto N = std::make_unique<lagrange_matrix_base_t[]>(n * n);
  auto D = std::make_unique<lagrange_matrix_base_t[]>(n * n);
  auto ipiv = std::make_unique<int[]>(n);

  // how to set t?
  double inf_norm =
      LAPACKE_dlange(CblasRowMajor, 'I', rows, rows, A.get(), leading_dim);
  assert(inf_norm > 0.0);
  int At_norm = static_cast<int>(inf_norm * t);
  int scale_exp = std::min(30, std::max(0, 1 + At_norm));

  double Ascal = t / std::pow(2.0, scale_exp);
  assert(Ascal > 0.0);
  cblas_dscal(n * n, Ascal, A.get(), 1);

  constexpr int q = 3;
  double c = 0.5;
  double sign = -1.0;

  cblas_dcopy(n * n, A.get(), 1, X_1.get(), 1);

  for (size_t i = 0; i < n * n; i++) {
    X_2[i] = 0.0;
    N[i] = 0.0;
    D[i] = 0.0;
  }

  for (size_t i = 0; i < static_cast<size_t>(rows); i++) {
    MIDX(N, i, i, n) = 1.0;
    MIDX(D, i, i, n) = 1.0;
  }

  cblas_daxpy(n * n, c, X_1.get(), 1, N.get(), 1);
  cblas_daxpy(n * n, sign * c, X_1.get(), 1, D.get(), 1);

  for (int i = 2; i <= q;) {
    c = c * (q - i + 1) / (i * (2 * q - i + 1));
    sign *= -1.0;

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rows, rows, rows,
                1.0, A.get(), leading_dim, X_1.get(), leading_dim, 0.0,
                X_2.get(), leading_dim);
    cblas_daxpy(n * n, c, X_2.get(), 1, N.get(), 1);
    cblas_daxpy(n * n, sign * c, X_2.get(), 1, D.get(), 1);

    i += 1;

    if (i > q) { break; }

    c = c * (q - i + 1) / (i * (2 * q - i + 1));
    sign *= -1.0;

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rows, rows, rows,
                1.0, A.get(), leading_dim, X_2.get(), leading_dim, 0.0,
                X_1.get(), leading_dim);
    cblas_daxpy(n * n, c, X_1.get(), 1, N.get(), 1);
    cblas_daxpy(n * n, sign * c, X_1.get(), 1, D.get(), 1);

    i += 1;
  }

  ret = LAPACKE_dgesv(CblasRowMajor, rows, rows, D.get(), leading_dim,
                      ipiv.get(), N.get(), leading_dim);
  assert(ret == 0);

  auto r1 = std::move(N);
  auto r2 = std::move(D);
  // N and D are invalid

  for (int i = 0; i < scale_exp; ++i) {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rows, rows, rows,
                1.0, r1.get(), leading_dim, r1.get(), leading_dim, 0.0,
                r2.get(), leading_dim);
    std::swap(r1, r2);
  }

  if (transposed) {
    cblas_domatcopy(CblasRowMajor, CblasTrans, rows, rows, 1.0, r1.get(),
                    leading_dim, r2.get(), leading_dim);
    for (int i = 0; i < rows; i++) { r2.get()[i * n] = 0.0; }
    r2.get()[0] = 1.0;
    std::swap(r1, r2);
  }

  return r1;
}

// first calculate Krylov subspace, but then calculate the expm
// operation of lower dimension with Pade
void expm_multiply_arnoldi_pade(const std::shared_ptr<Workspace> ws,
                                size_t rate_matrix, size_t clv_src,
                                size_t clv_dst, bool transposed,
                                size_t prob_matrix, double t) {
  //
  // allocate buffers
  //

  const int m = 20;
  const int rows = static_cast<int>(ws->matrix_rows());
  const int n = rows;
  int leading_dim = static_cast<int>(ws->leading_dimension());

  assert(rows > 0);
  assert(leading_dim > 0);

  auto A = std::make_unique<lagrange_matrix_base_t[]>(ws->matrix_size());

  auto e1_c = std::make_unique<lagrange_complex_t[]>(m);
  auto Hm_minus_theta_I_c =
      std::make_unique<lagrange_complex_t[]>(m * m * N_ROOT);
  auto y_c = std::make_unique<lagrange_complex_t[]>(m * N_ROOT);

  auto eHme1 = std::make_unique<lagrange_matrix_base_t[]>(m);
  auto H = std::make_unique<lagrange_matrix_base_t[]>((m + 1) * m);
  auto Q = std::make_unique<lagrange_matrix_base_t[]>(n * (m + 1));
  auto e1 = std::make_unique<lagrange_matrix_base_t[]>(m);

  auto ipiv_m = std::make_unique<int[]>(n);

  //
  // initialization
  //

  for (size_t i = 0; i < ws->matrix_size(); i++) {
    A.get()[i] = ws->rate_matrix(rate_matrix)[i] * t;
  }

  if (transposed) {
    cblas_dimatcopy(CblasRowMajor, CblasTrans, rows, rows, 1.0, A.get(),
                    leading_dim, leading_dim);
  }

  std::fill(e1_c.get(), e1_c.get() + m, 0.0);
  e1_c[0] = 1.0;

  std::fill(e1.get(), e1.get() + m, 0.0);
  e1[0] = 1.0;

  //
  // arnoldi iteration
  //

  std::fill(H.get(), H.get() + (m + 1) * m, 0.0);
  std::fill(Q.get(), Q.get() + n * (m + 1), 0.0);

  const double beta = cblas_dnrm2(n, ws->clv(clv_src), 1);

  cblas_dcopy(n, ws->clv(clv_src), 1, Q.get(), m + 1);
  cblas_dscal(n, 1.0 / beta, Q.get(), m + 1);

  for (int j = 0; j < m; j++) {
    // candidate for the next basis: w = A * v_j+1
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, A.get(), n, Q.get() + j,
                m + 1, 0.0, Q.get() + j + 1, m + 1);

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

  // pade scaling and squaring
  auto eHm = expm_pade(H, m, m, false, 1.0);

  cblas_dgemv(CblasRowMajor, CblasNoTrans, m, m, 1.0, eHm.get(), m, e1.get(), 1,
              0.0, eHme1.get(), 1);

  // update target clv directly
  cblas_dgemv(CblasRowMajor, CblasNoTrans, n, m, beta, Q.get(), m + 1,
              eHme1.get(), 1, 0.0, ws->clv(clv_dst), 1);
}

void expm_multiply_pade(const std::shared_ptr<Workspace> ws,
                        size_t rate_matrix_index, size_t clv_src,
                        size_t clv_dst, bool transposed,
                        size_t prob_matrix_index, double _t) {
  auto A = std::make_unique<lagrange_matrix_base_t[]>(ws->matrix_size());

  int ret = 0;
  int rows = static_cast<int>(ws->matrix_rows());
  int leading_dim = static_cast<int>(ws->leading_dimension());

  assert(rows > 0);
  assert(leading_dim > 0);

  for (size_t i = 0; i < ws->matrix_size(); i++) {
    A.get()[i] = ws->rate_matrix(rate_matrix_index)[i];
  }

  auto eA = expm_pade(A, rows, leading_dim, transposed, _t);

  // ws->update_prob_matrix(prob_matrix_index, eA.get());
  // ws->advance_clock();

  cblas_dgemv(CblasRowMajor, CblasNoTrans, ws->restricted_state_count(),
              ws->restricted_state_count(), 1.0, eA.get(),
              ws->leading_dimension(), ws->clv(clv_src), 1, 0.0,
              ws->clv(clv_dst), 1);
}

}  // namespace arnoldi

#endif  // LAGRANGE_CPP_ARNOLDI_H
