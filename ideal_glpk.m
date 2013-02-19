function ideal = ideal_glpk(constraints, lpr)
% IDEAL_GLPK  returns the ideal point of the maximal dominated coherent lower previsions
%
% Synopsis:
%    ideal = ideal_glpk(constraints, lpr)
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
%    ideal = a column vector with lenth the number of rows of constraints.A
%
% Background & Method:
%    Given a set of constraints [A, b] that describe the polytope of
%    coherent lower previsions, we can calculate the ideal point F of the
%    set of maximal coherent lower previsions dominated by a given lower
%    prevision P using a linear program for each component f (column of K):
%
%       maximize F(f) subject to A[F; μ] <= b and F <= P,
%
%    where μ may be a vector of additional variables in the constraint
%    definition; it is asumed that the first columns of A pertain to F.
%
% See also GLPK.

  m = length(lpr);
  [l, k] = size(constraints.A);

  % preallocate
  ideal = lpr;

  % prepare the common parts of the linear programs
  C = [eye(m), zeros(m, k - m)];
  A = [constraints.A; eye(m), zeros(m, k - m)];
  b = [constraints.B; lpr];
  ctype = repmat('U', 1, l + m);

  % for each gamble, complete the preparations for and solve the linear
  % program
  for k = 1:m
    c = C(k,:);
    [~, ideal(k), ~, ~] = glpk(c, A, b, [], [], ctype, [], -1);
  end

end