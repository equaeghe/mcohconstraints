function constraints = coh_free_constraints(K)
% COH_FREE_CONSTRAINTS  returns the 'free' constraints for coherence
%
% Synopsis:
%    constraints = coh_free_constraints(K)
%
% Input:
%    K = a nonnegative matrix with nonconstant columns ("gambles")
%
% Output:
%    constraints = a struct describing constraints P'λ <= α with three
%                  fields:
%                    * A, a matrix containing the constraint coefficients
%                         as rows
%                    * B, a column vector containing the corresponding
%                         constraint constants as components
%                    * lin, a column vector with indices of constraints
%                           that are actually equalities (== instead of <=)
%
% Background & Method:
%    To be coherent, a lower prevision defined on a set of gambles K must
%    belong to the polyhedron defined by the constraints SP <= SK'μ_S,
%    μ_S >= 0, and 1'μ_S == 1 for all matrices S that differ from the
%    identity matrix by at most one signchange.
%
% See also COH_FREE_CONSTRAINTS_FILE, ASL_FREE_CONSTRAINTS

  [n, m] = size(K);
  N = 1:n;
  M = 1:m;

  % preallocate
  A = zeros((m + 1) * (m + n + 2), m * (1 + n));
  b = zeros((m + 1) * (m + n + 2), 1);

  A(M,M) = eye(m);
  A(M,m+N) = -K';
  A((m+1)*m+N,m+N) = -eye(n);
  A((m+1)*(m+n)+1,m+N) = ones(1, n);
  b((m+1)*(m+n)+1) = 1;
  A((m+1)*(m+n+1)+1,m+N) = -ones(1, n);
  b((m+1)*(m+n+1)+1) = -1;

  for k = M
    S = eye(m);
    S(k,k) = -1;

    A(k*m+M,M) = S;
    A(k*m+M,m+k*n+N) = -S*K';
    A((m+1)*m+k*n+N,m+k*n+N) = -eye(n);
    A((m+1)*(m+n)+k+1,m+k*n+N) = ones(1, n);
    b((m+1)*(m+n)+k+1) = 1;
    A((m+1)*(m+n+1)+k+1,m+k*n+N) = -ones(1, n);
    b((m+1)*(m+n+1)+k+1) = -1;
  end

  constraints = struct('lin', [], 'A', A, 'B', b);

end