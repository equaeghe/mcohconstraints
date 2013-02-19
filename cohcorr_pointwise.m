function cohcorr = cohcorr_pointwise(constraints, lpr)
% COHCORR_POINTWISE  returns an upper approximation of the nadir point
%
% Synopsis:
%    cohcorr = cohcorr_pointwise(constraints, lpr)
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
%
% Background & Method:
%    EXPERIMENTAL STUFF, MAY BE BROKEN OR PRODUCe RESULTS THAT MAKE NO
%    SENSE WHATSOEVER.
%
% See also GLPK, COHCORR_BENSOLVE,
%          ASL_CONSTRAINTS_CDDMEX, COH_CONSTRAINTS_BENSOLVE.

  eps = 10^(-10);

  m = length(lpr);
  [l, k] = size(constraints.A);

  % use the ideal point as a starting point
  ideal = ideal_glpk(constraints, lpr);
  % preallocate
  cohcorr = ideal;

  % we only correct components whose ideal value is already lower,
  % THIS IS STUPID
  criterion = abs(lpr - ideal) > eps;
  indices = 1:m;
  indices = indices(criterion);

  % prepare the common parts of the linear programs
  C = [eye(m), zeros(m, k - m)];
  A = [constraints.A; eye(m), zeros(m, k - m)];
  b = [constraints.B; ideal];

  % for each gamble, complete the preparations for and solve the linear
  % program
  for n = indices
    c = C(n,:);
    ctype = [repmat('U', 1, l), repmat('S', 1, m)];
    ctype(l+n) = 'U';
    [~, cohcorr(n), ~, ~] = glpk(c, A, b, [], [], ctype, [], -1);
  end

end