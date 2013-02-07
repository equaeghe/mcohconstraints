function constraints = asl_constraints_bensolve(K)
% ASL_CONSTRAINTS_BENSOLVE  returns avoiding sure loss constraints
%
% Synopsis:
%    constraints = asl_constraints_bensolve(K)
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
%    So, to find the set of λ that characterize avoiding sure loss, we need
%    to solve the MOLP problem to the extreme points of the efficient set
%    and add positivity constraints that express our assumption P >= 0. The
%    MOLP solver used is bensolve.
%
% See also BENSOLVE, ASL_CONSTRAINTS_CDDMEX, COH_CONSTRAINTS_BENSOLVE.

  [n, m] = size(K);

  % solve the MOLP max Kλ s.t. Kλ <= 1 and λ >= 0,
  % i.e., min -λ s.t. [-K; I)λ >= [-1; 0] because our solver only returns
  % optimal objective vectors, not optimization vectors; we can make this
  % change because K is nonnegative:
  [~, ~, ~, lambdas, ~, ~] = bensolve(-eye(m), [-K; eye(m)], ...
                                      [-ones(n, 1); zeros(m, 1)]);

  % put into output form and add positivity constraints
  constraints.A = [-lambdas'; -eye(m)];
  constraints.B = [ones(size(lambdas, 2), 1); zeros(m, 1)];
  constraints.lin = [];

end