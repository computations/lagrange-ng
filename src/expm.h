/*!
 *   \file expm.hpp
 *
 *   Implement matrix exponential using pade approximation.
 *
 *  Copyright (c) 2007
 *  \author Tsai, Dung-Bang
 *
 */
//  Department of Physics,
//  National Taiwan University.
//
//  E-Mail : dbtsai [_at_] dbtsai [_dot_] org
//  Begine : 2007/11/20
//  Last modify : 2007/11/26
//  Version : v0.4
//
//  expm_pad computes the matrix exponential exp(H) for general matrixs,
//  including complex and real matrixs using the irreducible (p,p) degree
//  rational Pade approximation to the exponential
//  exp(z) = r(z) = (+/-)( I+2*(Q(z)/P(z))).
//
//  Usage :
//
//   U = expm_pad(H)
//   U = expm_pad(H, t),
//   U = expm_pad(H, t, p),
//
//  where t is a real number which is default set to 1.0 such that U=exp(t*H),
//  and p is internally set to 6 (recommended and gererally satisfactory).
//
//  See also MATLAB supplied functions, EXPM and EXPM1.
//
//  Reference :
//  EXPOKIT, Software Package for Computing Matrix Exponentials.
//  ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
//
// Use, modification and distribution are subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt).
//

#ifndef _BOOST_UBLAS_EXPM_
#define _BOOST_UBLAS_EXPM_
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/traits.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace boost {
namespace numeric {
namespace ublas {

template <typename MATRIX>
MATRIX
expm_pad(const MATRIX &H,
         typename type_traits<typename MATRIX::value_type>::real_type t = 1.0,
         const int p = 6) {
  typedef typename MATRIX::value_type value_type;
  typedef typename MATRIX::size_type size_type;
  typedef typename type_traits<value_type>::real_type real_value_type;
  assert(H.size1() == H.size2());
  assert(p >= 1);
  const size_type n = H.size1();
  const identity_matrix<value_type> I(n);
  matrix<value_type> U(n, n), H2(n, n), P(n, n), Q(n, n);
  real_value_type norm = 0.0;
  // Calcuate Pade coefficients
  vector<real_value_type> c(p + 1);
  c(0) = 1;
  for (size_type i = 0; i < p; ++i)
    c(i + 1) = c(i) * ((p - i) / ((i + 1.0) * (2.0 * p - i)));
  // Calcuate the infinty norm of H, which is defined as the largest row sum of
  // a matrix
  for (size_type i = 0; i < n; ++i) {
    real_value_type temp = 0.0;
    for (size_type j = 0; j < n; j++)
      temp += std::abs(H(i, j));
    norm = t * std::max<real_value_type>(norm, temp);
  }
  // If norm = 0, and all H elements are not NaN or infinity but zero,
  // then U should be identity.
  if (norm == 0.0) {
    bool all_H_are_zero = true;
    for (size_type i = 0; i < n; i++)
      for (size_type j = 0; j < n; j++)
        if (H(i, j) != value_type(0.0))
          all_H_are_zero = false;
    if (all_H_are_zero == true)
      return I;
    // Some error happens, H has elements which are NaN or infinity.
    std::cerr << "Null input error in the template expm_pad.\n";
    // std::cout << "Null INPUT : " << H <<"\n";
    exit(0);
  }
  // Scaling, seek s such that || H*2^(-s) || < 1/2, and set scale = 2^(-s)
  int s = 0;
  real_value_type scale = 1.0;
  if (norm > 0.5) {
    s = std::max<int>(0, static_cast<int>((log(norm) / log(2.0) + 2.0)));
    scale /= real_value_type(std::pow(2.0, s));
    U.assign((scale * t) *
             H); // Here U is used as temp value due to that H is const
  } else
    U.assign(H);

  // Horner evaluation of the irreducible fraction, see the following ref above.
  // Initialise P (numerator) and Q (denominator)
  H2.assign(prod(U, U));
  Q.assign(c(p) * I);
  P.assign(c(p - 1) * I);
  size_type odd = 1;
  for (size_type k = p - 1; k > 0; --k) {
    (odd == 1) ? (Q = (prod(Q, H2) + c(k - 1) * I))
               : (P = (prod(P, H2) + c(k - 1) * I));
    odd = 1 - odd;
  }
  (odd == 1) ? (Q = prod(Q, U)) : (P = prod(P, U));
  Q -= P;
  // In origine expokit package, they use lapack ZGESV to obtain inverse matrix,
  // and in that ZGESV routine, it uses LU decomposition for obtaing inverse
  // matrix. Since in ublas, there is no matrix inversion template, I simply use
  // the build-in LU decompostion package in ublas, and back substitute by
  // myself.

  // Implement Matrix Inversion
  permutation_matrix<size_type> pm(n);
  int res = lu_factorize(Q, pm);
  if (res != 0) {
    std::cerr << "Matrix inversion error in the template expm_pad.\n";
    exit(0);
  }
  // H2 is not needed anymore, so it is temporary used as identity matrix for
  // substituting.
  H2.assign(I);
  lu_substitute(Q, pm, H2);
  (odd == 1) ? (U.assign(-(I + real_value_type(2.0) * prod(H2, P))))
             : (U.assign(I + real_value_type(2.0) * prod(H2, P)));
  // Squaring
  for (size_type i = 0; i < s; ++i)
    U = (prod(U, U));
  return U;
}
} // namespace ublas
} // namespace numeric
} // namespace boost
#endif
