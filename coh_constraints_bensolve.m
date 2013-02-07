function constraints = coh_constraints_bensolve(K)
% COH_CONSTRAINTS_BENSOLVE  returns coherence constraints
%
% Synopsis:
%    constraints = coh_constraints_bensolve(K)
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
%    So, to find the set of λ that characterize avoiding sure loss, we need
%    to the MOLP problem to the extreme points of the efficient set for
%    both α = 0 and α = 1, and add positivity constraints that express our
%    assumption P >= 0.
%
%    However, λ ~>= 0 is not a linear constraint, so we must work case by
%    case: For each possibly negative component λ_k, we solve the MOLP
%    max Kλ s.t. λ >= 0 except λ_k <= 0 and Kλ <= α. (We must also not
%    forget the case λ >= 0 and Kλ <= 1.) To finish, we need to group the
%    constraints and remove any redundant ones. The MOLP solver used is
%    bensolve; the redundancy remover is cddlib through cddmex.
%
% See also BENSOLVE, CDDMEX, COH_CONSTRAINTS_CDDMEX,
%          ASL_CONSTRAINTS_BENSOLVE.

  [n, m] = size(K);

  % preallocate some datastructures
  lambdas0 = cell(1,m);
  lambdas1 = cell(1,m+1);
  
  % case λ >= 0 except for one component λ_k <= 0:
  for k = 1:m
    % create diagonal matrix indicating nonpositive component λ_k:
    S = eye(m);
    S(k,k) = -1;

    % case α = 0:
    % solve the MOLP max Kλ s.t. Kλ <= 0 and 0 <= λ <= 1 except
    % -1 <= λ_k <= 0, i.e., min -λ s.t. [-K; S; -S]λ >= [0; 0; -1] because
    % our solver only returns optimal objective vectors, not optimization
    % vectors; we can make this change because K is nonnegative:
    [~, ~, ~, lambdas0{k}, ~, ~] = bensolve(-eye(m), [-K; S; -S], ...
                                            [zeros(n, 1); zeros(m, 1); ...
                                             -ones(m, 1)]);
    lambdas0{k} = -lambdas0{k}';

    % case α = 1:
    % solve the MOLP max Kλ s.t. Kλ <= 1 and λ >= 0 except λ_k <= 0,
    % i.e., min -λ s.t. [-K; S]λ >= [-1; 0] because our solver only returns
    % optimal objective vectors, not optimization vectors; we can make this
    % change because K is nonnegative:
    [~, ~, ~, lambdas1{k}, ~, ~] = bensolve(-eye(m), [-K; S], ...
                                            [-ones(n, 1); zeros(m, 1)]);
    lambdas1{k} = -lambdas1{k}';
  end

  % case λ >= 0 and α = 1:
  % solve the MOLP max Kλ s.t. Kλ <= 1 and λ >= 0,
  % i.e., min -λ s.t. [-K; I)λ >= [-1; 0] because our solver only returns
  % optimal objective vectors, not optimization vectors; we can make this
  % change because K is nonnegative:
  [~, ~, ~, lambdas1{m+1}, ~, ~] = bensolve(-eye(m), [-K; eye(m)], ...
                                            [-ones(n, 1); zeros(m, 1)]);
  lambdas1{m+1} = -lambdas1{m+1}';

  % group constraints of the same type (α = 0 or α = 1) and remove
  % recurring ones
  lambdas0 = unique(vertcat(lambdas0{1:end}), 'rows');
  lambdas1 = unique(vertcat(lambdas1{1:end}), 'rows');

  % add positivity constraints and then remove redundant constraints
  constraints = cddmex('reduce_h', ...
                       struct('A', [lambdas0; lambdas1; -eye(m)], ...
                              'B', [zeros(size(lambdas0, 1), 1); ...
                                    ones(size(lambdas1, 1), 1); ...
                                    zeros(m, 1)]));

end