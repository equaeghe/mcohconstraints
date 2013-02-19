function natex = natex_bensolve(constraints, lpr)
% NATEX_BENSOLVE  returns the natural extension using a MOLP solver
%
% Synopsis:
%    natex = natex_bensolve(constraints, lpr)
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
%    natex = a column vector with lenth the number of rows of constraints.A
%
% Background & Method:
%    Given a set of constraints [A, b] that describe the polytope of
%    coherent lower previsions, we can calculate the natural extension E of
%    any lower prevision P using the MOLP
%
%       minimize E subject to A[E; μ] <= b and E >= P,
%
%    where μ may be a vector of additional variables in the constraint
%    definition; it is asumed that the first columns of A pertain to E.
%
% See also COH_CONSTRAINTS_CDDMEX, COH_FREE_CONSTRAINTS, BENSOLVE,
%          NATEX_DIRECT, NATEX_LENV.

  m = length(lpr);
  k = size(constraints.A, 2);

  % prepare the MOLP
  C = [eye(m), zeros(m, k-m)];
  A = [-constraints.A; eye(m), zeros(m, k-m)];
  b = [-constraints.B; lpr];
  options = struct('info', 0);

  % solve the MOLP
  [~, ~, ~, natex, ~, ~] = bensolve(C, A, b, [], [], [], options);

end