function natex = natex_direct(K, lpr)
% NATEX_DIRECT  returns the natural extension using the direct LP approach
%
% Synopsis:
%    natex = natex_direct(K, lpr)
%
% Input:
%    K = a matrix ("gambles")
%    lpr = a column vector with length the number of rows of K
%
% Output:
%    natex = a column vector with lenth the number of rows of K
%
% Background & Method:
%    For each gamble f (column of K), the natural extension E(f) of a
%    lower prevision P can be calculated using the following linear
%    program:
%
%       E(f) = max α subject to (K-P)λ + α <= f and λ >= 0
%
% See also GLPK, NATEX_BENSOLVE, NATEX_LENV.

  [n, m] = size(K);
  natex = lpr;

  % prepare the common parts of the linear programs
  c = [zeros(m, 1); 1; -1];
  A = [K - repmat(lpr', n, 1), ones(n, 1), -ones(n, 1)];
  ctype = repmat('U', 1, n);

  % for each gamble, complete the preparations for and solve the linear
  % program
  for k = 1:m
    b = K(:, k);
    [~, natex(k), ~, ~] = glpk(c, A, b, [], [], ctype, [], -1);
  end

end