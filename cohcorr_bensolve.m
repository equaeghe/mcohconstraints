function [cohcorr, maxnum] = cohcorr_bensolve(constraints, lpr)

  m = length(lpr);
  k = size(constraints.A, 2);

  C = [-eye(m), zeros(m, k-m)];
  A = [-constraints.A; -eye(m), zeros(m, k-m)];
  b = [-constraints.B; -lpr];
  options = struct('info', 0, 'vert_enum', 'C');

  [~, ~, ~, cohcorr_maximals, ~, ~] = bensolve(C, A, b, [], [], [], options);

  cohcorr = min(-cohcorr_maximals')';
  maxnum = size(cohcorr_maximals, 2);

end