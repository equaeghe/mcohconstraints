function natex = natex_bensolve(constraints, lpr)

  m = length(lpr);
  k = size(constraints.A, 2);

  C = [eye(m), zeros(m, k-m)];
  A = [-constraints.A; eye(m), zeros(m, k-m)];
  b = [-constraints.B; lpr];
  options = struct('info', 0);

  [~, ~, ~, natex, ~, ~] = bensolve(C, A, b, [], [], [], options);

end