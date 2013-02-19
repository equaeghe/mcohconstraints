function [cohcorr, maxnum] = cohcorr_bensolve(constraints, lpr)
% COHCORR_BENSOLVE  returns the nadir point of of the maximal dominated coherent lower previsions
%
% Synopsis:
%    [cohcorr, maxnum] = cohcorr_bensolve(constraints, lpr)
%
% Input:
%    constraints = a struct describing coherence constraints P'λ <= α with
%                  three fields (it may be there are more optimization
%                  vector components than just P):
%                    * A, a matrix containing the constraint coefficients λ
%                         as rows
%                    * B, a column vector containing the corresponding
%                         constraint constants α as components
%                    * lin, a column vector with indices of constraints
%                           that are actually equalities (== instead of <=)
%    lpr = a column vector with lenth the number of rows of constraints.A
%
% Output:
%    cohcorr = a column vector with lenth the number of rows of
%              constraints.A
%    maxnum = the number of extreme maximal dominated coherent lower
%             previsions
%
% Background & Method:
%    Given a set of constraints [A, b] that describe the polytope of
%    coherent lower previsions, we can calculate the extreme maximal
%    coherent lower previsions Q dominated by P using the MOLP
%
%       maximize Q subject to A[Q; μ] <= b and Q <= P,
%
%    where μ may be a vector of additional variables in the constraint
%    definition; it is asumed that the first columns of A pertain to Q.
%    Then the nadir point D is the lower envelope of all the extreme
%    maximal solutions of the MOLP.
%
% See also BENSOLVE, ASL_CONSTRAINTS_CDDMEX, COH_CONSTRAINTS_BENSOLVE.

  m = length(lpr);
  k = size(constraints.A, 2);

  % prepare the MOLP
  C = [-eye(m), zeros(m, k-m)];
  A = [-constraints.A; -eye(m), zeros(m, k-m)];
  b = [-constraints.B; -lpr];
  options = struct('info', 0, 'vert_enum', 'C');

  % solve the MOLP
  [~, ~, ~, cohcorr_maximals, ~, ~] = bensolve(C, A, b, [], [], [], ...
                                               options);

  % take the lower envelope (sign compensation due to bensolve sign
  % conventions)
  cohcorr = min(-cohcorr_maximals, [], 2);
  maxnum = size(cohcorr_maximals, 2);

end