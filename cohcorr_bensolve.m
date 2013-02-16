function [cohcorr, maxnum] = cohcorr_bensolve(constraints, lpr)

  m = length(lpr);

  [~, ~, ~, cohcorr_maximals, ~, ~] = bensolve(-eye(m), ...
                                               [-constraints.A; -eye(m)], ...
                                               [-constraints.B; -lpr], ...
                                               [], [], [], struct('info', 0));

  cohcorr = min(-cohcorr_maximals')';
  maxnum = size(cohcorr_maximals, 2);

end