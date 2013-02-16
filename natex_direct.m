function natex = natex_direct(K, lpr)

  [n, m] = size(K);
  natex = lpr;

  c = [zeros(m, 1); 1; -1];
  A = [K - repmat(lpr', n, 1), ones(n, 1), -ones(n, 1)];

  for k = 1:m
    b = K(:, k);
    [~, natex(k), ~, ~] = glpk(c, A, b, [], [], repmat('U', 1, n), [], -1);
  end

end