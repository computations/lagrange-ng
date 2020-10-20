#ifndef _LAGRANGE_OPERATION_H
#define _LAGRANGE_OPERATION_H

#include "Workspace.h"
#include <blaze/Math.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <memory>

class ExpmOperation {
public:
  ExpmOperation(size_t prob_matrix, size_t rate_matrix, double t)
      : _prob_matrix_index{prob_matrix},
        _rate_matrix_index{rate_matrix}, _t{t} {}

  void eval(std::shared_ptr<Workspace> ws);

private:
  size_t _prob_matrix_index;
  size_t _rate_matrix_index;
  double _t;
};

class DispersionOperation {
public:
  DispersionOperation(size_t top, size_t bot, double brlen, size_t prob_matrix,
                      size_t rate_matrix)
      : _top_clv{top}, _bot_clv{bot}, _prob_matrix_index{prob_matrix},
        _expm_op{new ExpmOperation{prob_matrix, rate_matrix, brlen}} {}

  DispersionOperation(size_t top, size_t bot, size_t prob_matrix,
                      std::shared_ptr<ExpmOperation> op)
      : _top_clv{top}, _bot_clv{bot},
        _prob_matrix_index{prob_matrix}, _expm_op{op} {}

  DispersionOperation(size_t top, size_t bot, size_t prob_matrix)
      : DispersionOperation(top, bot, prob_matrix, nullptr) {}

  void eval(std::shared_ptr<Workspace> ws);

private:
  /* Remember, the top and bottom clv indexes could be the same. This is to save
   * on storage when computing different periods along a single branch. The idea
   * is that we just apply the matrix multiplication over and over again to the
   * same CLV. Matvecs can't be done in place, so there will need to be
   * temporary storage if we do that.
   */
  size_t _top_clv;
  size_t _bot_clv;
  size_t _prob_matrix_index;

  /* _expm_op is an owning pointer, and can be null */
  std::shared_ptr<ExpmOperation> _expm_op;
};

class SplitOperation {
public:
  SplitOperation(size_t lchild_clv, size_t rchild_clv, size_t parent_clv,
                 DispersionOperation *lops, size_t lops_count,
                 DispersionOperation *rops, size_t rops_count)
      : _lbranch_clv_index{lchild_clv}, _rbranch_clv_index{rchild_clv},
        _parent_clv_index{parent_clv}, _lbranch_ops{lops},
        _lbranch_ops_count{lops_count}, _rbranch_ops{rops}, _rbranch_ops_count{
                                                                rops_count} {}

  SplitOperation(size_t lchild_clv_top, size_t rchild_clv_top,
                 size_t lchild_clv_bot, size_t rchild_clv_bot, double lbrlen,
                 double rbrlen, size_t lprob_mat, size_t rprob_mat,
                 size_t lrate_matrix, size_t rrate_matrix, size_t parent_clv)
      : _lbranch_clv_index{lchild_clv_top}, _rbranch_clv_index{rchild_clv_top},
        _parent_clv_index{parent_clv}, _lbranch_ops{new DispersionOperation{
                                           lchild_clv_top, lchild_clv_bot,
                                           lbrlen, lprob_mat, lrate_matrix}},
        _lbranch_ops_count{1}, _rbranch_ops{new DispersionOperation{
                                   rchild_clv_top, rchild_clv_bot, rbrlen,
                                   rprob_mat, rrate_matrix}},
        _rbranch_ops_count{1} {}

  SplitOperation(size_t lchild_clv_top, size_t lchild_clv_bot,
                 size_t rchild_clv_top, size_t rchild_clv_bot, double lbrlen,
                 double rbrlen, size_t prob_mat, size_t rate_matrix,
                 size_t parent_clv)
      : SplitOperation(lchild_clv_top, rchild_clv_top, lchild_clv_bot,
                       rchild_clv_bot, lbrlen, rbrlen, prob_mat, prob_mat,
                       rate_matrix, rate_matrix, parent_clv) {}

  void eval(std::shared_ptr<Workspace>);

private:
  std::vector<lagrange_region_split_t> generate_splits(uint64_t state,
                                                       size_t regions);

  size_t _lbranch_clv_index;
  size_t _rbranch_clv_index;
  size_t _parent_clv_index;

  std::shared_ptr<DispersionOperation> _lbranch_ops;
  size_t _lbranch_ops_count;
  std::shared_ptr<DispersionOperation> _rbranch_ops;
  size_t _rbranch_ops_count;
};

#endif
