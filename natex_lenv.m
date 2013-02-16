function natex = natex_lenv(K, lpr)

  [n, m] = size(K);
  natex = lpr;

  A = [K'; ones(1, n)];
  b = [lpr; 1];

  for k = 1:m
    c = K(:, k);
    [~, natex(k), ~, ~] = glpk(c, A, b, [], [], [repmat('L', 1, m), 'S']);
  end

end