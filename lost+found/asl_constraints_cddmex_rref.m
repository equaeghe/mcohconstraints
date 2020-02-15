function [lambdas, alphas] = asl_constraints_cddmex_rref(K)
% asl_constraints_cddmex_rref  returns a (matrix,vector)-pair of 
%                              coefficients λ and constant α defines a
%                              constraint P'λ <= α that a lower prevision P
%                              needs to satisfy to avoid sure loss
%
% Synopsis:
%    [lambdas, alphas] = asl_constraints_cddmex_rref(K)
%
% Input:
%    K = a matrix whose columns are interpreted as gambles
%
% Output:
%    lambdas = a matrix whose columns are constraint coefficients
%    alphas = a matrix whose components are constraint constants
%
% Remark:
%    This function makes use of Gauss-Jordan elimination under the hood.
%    For some matrices, this may introduce more round-off errors than other
%    methods would.

% Method:
%    To avoid sure loss, a lower prevision P defined on a set of gambles K
%    must (pointwise) satisfy P'λ <= max(Kλ) for all vectors λ >= 0. More
%    stringent constraints correspond to (pointwise) higher values of λ. If
%    we add slack variables σ >= 0, we can enforce max(Kλ) = Kλ + Iσ = α,
%    for some constant α. So then the constraints become P'λ <= α.
%
%    So we get a constraint for every λ such that Kλ + Iσ = α for some σ.
%    The set of such λ is convex, so we can limit attention to the extreme
%    points of this convex set. To generate these, we could, crudely put,
%    perform vertex enumeration on {Kλ + Iσ = α, λ >= 0, σ >= 0}, project
%    onto λ-space and remove non-extreme projections.
%
%    However, because vertex enumeration is a computationally expensive
%    operation, we can first partially solve Kλ + Iσ = α. Let N be the
%    matrix with columns that form a basis of the (column/right) nullspace
%    of [K I], then we can write [λ; σ] = κ + Nμ, where κ is any one
%    solution of Kλ + Iσ = α and μ may be any real vector. Then, crudely
%    put, we only need to do vertex enumeration on {Nμ >= -κ}, whereby the
%    computational efficiency is improved because the size of μ is the same
%    as the size of λ, so it is smaller than the size of [λ; σ].
%
%    Normalization and nonnegativity of λ allows us to limit α to 0 and 1.
%    To get a minimal set of constraints, we need to remove possible
%    redundant constraints after combining the constraints for each case,
%    which are minimal on their own.

  [n, m] = size(K);
  % calculate a basis, indices to the basis vectors, and the nullspace of K
  % augmented with the identity matrix (for the slack variables)
  [B, b, ~, N] = matrixdigest_rref([K, eye(n)]);

  % case α = 0;
  % because Kλ + Iσ = 0 does not bound λ, we add bound constraints Nμ <= 1:
  kappa0 = zeros(m + n, 1);
  mus0 = cddmex('extreme', struct('A', [-N; N], ...
                                  'B', [kappa0; ones(m + n, 1)])).V;
  lambdas0 = mus0 * N(1:m,:)';
  clear kappa0 mus0;
  lambdas0 = cddmex('reduce_v', struct('V', unique(lambdas0, 'rows'))).V;

  % case α = 1:
  kappa1 = zeros(m + n, 1);
  kappa1(b) = B \ ones(n, 1);
  mus1 = cddmex('extreme', struct('A', -N, 'B', kappa1)).V;
  lambdas1 = repmat(kappa1(1:m)', size(mus1, 1), 1) + mus1 * N(1:m,:)';
  clear kappa1 mus1;
  lambdas1 = cddmex('reduce_v', struct('V', unique(lambdas1, 'rows'))).V;

  % add positivity constraints and then remove redundant constraints
  h = cddmex('reduce_h', struct('A', [lambdas0; lambdas1; -eye(m)], ...
                                'B', [zeros(size(lambdas0, 1), 1); ...
                                      ones(size(lambdas1, 1), 1); ...
                                      zeros(m, 1)]));
  clear lambdas0 lambdas1;

  lambdas = h.A';
  alphas = h.B';

end