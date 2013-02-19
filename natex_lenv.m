function natex = natex_lenv(K, lpr)
% NATEX_LENV  returns the natural extension with lower envelope LP approach
%
% Synopsis:
%    natex = natex_lenv(K, lpr)
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
%       E(f) = min f'μ subject to K'μ <= P and μ >= 0 and 1'μ = 1
%
% See also GLPK, NATEX_BENSOLVE, NATEX_DIRECT.

  [n, m] = size(K);
  natex = lpr;

  % prepare the common parts of the linear programs
  A = [K'; ones(1, n)];
  b = [lpr; 1];
  ctype = [repmat('L', 1, m), 'S'];

  % for each gamble, complete the preparations for and solve the linear
  % program
  for k = 1:m
    c = K(:, k);
    [~, natex(k), ~, ~] = glpk(c, A, b, [], [], ctype);
  end

end