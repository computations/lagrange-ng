/*
 * Created by sunflower on 24.07.22.
 * Last edited by Ben Bettisworth on 2022-08-30
 */

#ifndef LAGRANGE_CPP_ARNOLDI_H
#define LAGRANGE_CPP_ARNOLDI_H

#include <memory>

#include "Workspace.hpp"

namespace lagrange { namespace expm {

/*
 * Combination of Arnoldi and Chebyshev approximation
 * see Saad 1990 and Saad 1992
 */
void arnoldi_chebyshev(const std::shared_ptr<Workspace> ws,
                       size_t rate_matrix_index,
                       size_t clv_src_index,
                       size_t clv_dst_index,
                       bool transpose,
                       double t);

void arnoldi_pade(const std::shared_ptr<Workspace> ws,
                       size_t rate_matrix_index,
                       size_t clv_src_index,
                       size_t clv_dst_index,
                       bool transpose,
                       double t);
}}  // namespace lagrange::expm
#endif  // LAGRANGE_CPP_ARNOLDI_H
