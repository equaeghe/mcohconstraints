function constraints = asl_constraints_cddmex(K)
% ASL_CONSTRAINTS_CDDMEX  returns avoiding sure loss constraints
%
% Synopsis:
%    constraints = asl_constraints_cddmex(K)
%
% Input:
%    K = a nonnegative matrix with nonconstant columns ("gambles")
%
% Output:
%    constraints = a struct describing constraints P'λ <= α with three
%                  fields:
%                    * A, a matrix containing the constraint coefficients λ
%                         as rows
%                    * B, a column vector containing the corresponding
%                         constraint constants α as components
%                    * lin, a column vector with indices of constraints
%                           that are actually equalities (== instead of <=)
%
% Background & Method:
%    To avoid sure loss, a lower prevision P defined on a set of gambles K
%    must (pointwise) satisfy P'λ <= max(Kλ) for all vectors λ >= 0. Define
%    α = max(Kλ), then the constraints become P'λ <= α. More stringent
%    constraints correspond to (pointwise) higher values of λ if we assume
%    P nonnegative.
%
%    The MOLP problem pointwise-maximize Kλ subject to Kλ <= α and λ >= 0
%    has as solutions the maximally stringent λ; because K is nonnegative,
%    at least one component of Kλ attains the bound α. These solutions form
%    the so-called pareto efficient set, a connected union of convex sets.
%    The finite number of extreme points of the efficient set are 
%    sufficient to characterize avoiding sure loss, because the constraints
%    P'λ <= α are linear.
%
%    This MOLP problem can be solved by vertex enumeration: Observe that
%    the convex polytope {λ >= 0: Kλ <= α} is empty for α < 0 and reduces
%    to the zero vector for α == 0 because K is nonnegative with
%    nonconstant columns. So only the case α > 0 must be considered. We may
%    take α = 1, as λ can be rescaled. The extreme points of the efficient
%    set form a subset of the extreme points of {λ >= 0: Kλ <= 1}.
%
%    So, to find the set of λ that characterize avoiding sure loss, we need
%    to vertex enumerate {λ >= 0: Kλ <= 1}, add positivity constraints that
%    express our assumption P >= 0, and remove the redundant constraints
%    (non-efficient vertices) from the resulting set. The vertex enumerator
%    and redundancy remover is cddlib through cddmex.
%
% See also CDDMEX, ASL_CONSTRAINTS_BENSOLVE, COH_CONSTRAINTS_CDDMEX.

  [n, m] = size(K);

  % vertex enumerate {λ >= 0: Kλ <= 1}, i.e., {[K; -I)λ <= [1; 0]}:
  H = struct('A', [K; -eye(m)], 'B', [ones(n, 1); zeros(m, 1)]);
  V = getfield(cddmex('extreme', H), 'V');

  % add positivity constraints and then remove redundant constraints
  H = struct('A', [V; -eye(m)], ...
             'B', [ones(size(V, 1), 1); zeros(m, 1)]);
  constraints = cddmex('reduce_h', H);

end