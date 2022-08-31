/*
 * Created by sunflower on 24.07.22.
 * Last edited by Ben Bettisworth on 2022-08-30
 */

#ifndef LAGRANGE_CPP_ARNOLDI_H
#define LAGRANGE_CPP_ARNOLDI_H

#include <algorithm>
#include <memory>

#include "Common.h"
#include "Workspace.h"

/*
 * expm_multiply_arnoldi_chebyshev:
 */

namespace expm {
/*
 * Combination of Arnoldi and Chebyshev approximation
 * see Saad 1990 and Saad 1992
 */
void multiply_arnoldi_chebyshev(const std::shared_ptr<Workspace> ws,
                                size_t rate_matrix, size_t clv_src,
                                size_t clv_dst, bool transposed, double t);
}  // namespace expm

#endif  // LAGRANGE_CPP_ARNOLDI_H
