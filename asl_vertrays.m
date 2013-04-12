function constraints = asl_vertrays(K)
% ASL_VERTRAYS  returns avoiding sure loss constraints
%
% Synopsis:
%    asl_vertrays(K)
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
%    To avoid sure loss, a lower prevision defined on a set of gambles K
%    must be dominated by a linear prevision, i.e., a convex combination of
%    degenerate previsions; these correspond to the rows of K. This means
%    that it belongs to the polyhedron defined by the degenerate previsions
%    as vertices and negative directions as extreme rays.
%
%  See also CDDMEX, ASL_VERTRAYS_FILE, COH_VERTRAYS

  H = cddmex('hull', struct('V', K, 'R', -eye(length(K))));
  constraints = cddmex('reduce_h', H);

end