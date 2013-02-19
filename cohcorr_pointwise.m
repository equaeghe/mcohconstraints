function cohcorr = cohcorr_pointwise(constraints, lpr)

  eps = 10^(-10);

  m = length(lpr);
  [l, k] = size(constraints.A);

  ideal = ideal_glpk(constraints, lpr);
  cohcorr = ideal;

  criterion = abs(lpr - ideal) > eps;
  indices = 1:m;
  indices = indices(criterion);

  C = [eye(m), zeros(m, k - m)];
  A = [constraints.A; eye(m), zeros(m, k - m)];
  b = [constraints.B; ideal];

  for n = indices
    c = C(n,:);
    ctype = [repmat('U', 1, l), repmat('S', 1, m)];
    ctype(l+n) = 'U';
    [~, cohcorr(n), ~, ~] = glpk(c, A, b, [], [], ctype, [], -1);
  end

end