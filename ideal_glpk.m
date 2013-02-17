function ideal = ideal_glpk(constraints, lpr)

  m = length(lpr);
  [l, k] = size(constraints.A);
  ideal = lpr;

  C = [eye(m), zeros(m, k - m)];
  A = [constraints.A; eye(m), zeros(m, k - m)];
  b = [constraints.B; lpr];

  for k = 1:m
    c = C(k,:);
    [~, ideal(k), ~, ~] = glpk(c, A, b, [], [], repmat('U', 1, l + m), [], -1);
  end

end