function constraints = asl_free_constraints(K)
% ASL_FREE_CONSTRAINTS  returns the 'free' constraints for avoiding sure loss
%
% Synopsis:
%    constraints = asl_free_constraints(K)
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
%    To avoid sure loss, a lower prevision defined on a set of gambles K must
%    be dominated by a linear prevision, i.e., a convex combination of
%    degenerate previsions; these correspond to the rows of K. This means that
%    it belongs to the polyhedron defined by the constraints P <= K'μ, μ >= 0,
%    and 1'μ == 1.
%
% See also ASL_FREE_CONSTRAINTS_FILE, ASL_FREE_CONSTRAINTS

  [n, m] = size(K);

  A = [eye(m), -K'; zeros(n, m), -eye(n); ...
       zeros(1, m), ones(1, n); zeros(1, m), -ones(1, n)];
  b = [zeros(m + n, 1); 1; -1];

  constraints = struct('lin', [], 'A', A, 'B', b);

end