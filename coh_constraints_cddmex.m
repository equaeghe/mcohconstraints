function constraints = coh_constraints_cddmex(K)
% COH_CONSTRAINTS_CDDMEX  returns coherence constraints
%
% Synopsis:
%    constraints = coh_constraints_cddmex(K)
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
%    To be coherent, a lower prevision P defined on a set of gambles K must
%    (pointwise) satisfy P'λ <= max(Kλ) for all vectors λ with at most one 
%    nonpositive component; denote this λ ~>= 0. Define α = max(Kλ), then
%    the constraints become P'λ <= α. More stringent constraints correspond
%    to (pointwise) higher values of λ if we assume P nonnegative.
%
%    The MOLP problem pointwise-maximize Kλ subject to Kλ <= α and λ ~>= 0
%    has as solutions the maximally stringent λ; because K is nonnegative,
%    at least one component of Kλ attains the bound α. These solutions form
%    the so-called pareto efficient set, a connected union of convex sets.
%    The finite number of extreme points of the efficient set are
%    sufficient to characterize avoiding sure loss, because the constraints
%    P'λ <= α are linear.
%
%    This MOLP problem can be solved by vertex enumeration: Observe that
%    the convex polytope {λ >= 0: Kλ <= α} is empty for α < 0. So only the
%    case α >= 0 must be considered. We may take α = 0 with bounded λ or
%    α = 1, as λ can be rescaled. For each α, the extreme points of the
%    efficient set form a subset of the extreme points of
%    {λ >= 0: Kλ <= α}.
%
%    So, to find the set of λ that characterize avoiding sure loss, we need
%    to vertex enumerate {λ ~>= 0: Kλ <= α} for both α = 0 and α = 1, add
%    positivity constraints that express our assumption P >= 0, and remove
%    the redundant constraints (non-efficient vertices) from the resulting
%    set. The vertex enumerator and redundancy remover is cddlib through
%    cddmex.
%
%    However, λ ~>= 0 is not a linear constraint, so we must work case by
%    case: For each possibly negative component λ_k, we vertex enumerate
%    {λ >= 0 except λ_k <= 0: Kλ <= α}, or equivalently,
%    {λ >= 0: Kkλ <= α}, with Kk the matrix that is equal to K except with
%    negated column k. (We must also not forget the case
%    {λ >= 0: Kλ <= 1}.)
%
% See also CDDMEX, COH_CONSTRAINTS_BENSOLVE, ASL_CONSTRAINTS_CDDMEX.

  [n, m] = size(K);

  % preallocate some datastructures
  V0 = cell(1,m);
  V1 = cell(1,m+1);

  % case λ >= 0 except for one component λ_k <= 0:
  for k = 1:m
    % modify gamble matrix to simulate nonpositive λ_k
    Kk = K;
    Kk(:,k) = -K(:,k);

    % case α = 0:
    % vertex enumerate {0 <= λ <= 1: Kkλ <= 0},
    % i.e., {[Kk; -I; I)λ <= [0; 0; 1]}:
    H = struct('A', [Kk; -eye(m); eye(m)], ...
               'B', [zeros(n, 1); zeros(m, 1); ones(m, 1)]);
    V0{k} = getfield(cddmex('extreme', H), 'V');       

    % case α = 1:
    % vertex enumerate {λ >= 0: Kkλ <= 1}, i.e., {[Kk; -I)λ <= [1; 0]}:
    H = struct('A', [Kk; -eye(m)], 'B', [ones(n, 1); zeros(m, 1)]);
    V1{k} = getfield(cddmex('extreme', H), 'V');

    % modify lambda matrix to express nonpositive λ_k
    V0{k}(:,k) = -V0{k}(:,k);
    V1{k}(:,k) = -V1{k}(:,k);
  end

  % case λ >= 0 and α = 1:
  % vertex enumerate {λ >= 0: Kλ <= 1}, i.e., {[K; -I)λ <= [1; 0]}:
  H = struct('A', [K; -eye(m)], 'B', [ones(n, 1); zeros(m, 1)]);
  V1{m+1} = getfield(cddmex('extreme', H), 'V');

  % group constraints of the same type (α = 0 or α = 1) and remove
  % recurring ones
  V0 = unique(vertcat(V0{1:end}), 'rows');
  V1 = unique(vertcat(V1{1:end}), 'rows');

  % add positivity constraints and then remove redundant constraints
  H = struct('A', [V0; V1; -eye(m)], ...
             'B', [zeros(size(V0, 1), 1); ones(size(V1, 1), 1); ...
                   zeros(m, 1)]);
  constraints = cddmex('reduce_h', H);

end