function natex = natex_bensolve(constraints, lpr)

  m = length(lpr)

  [~, ~, ~, natex, ~, ~] = bensolve(eye(m), [-constraints.A; eye(m)], ...
                                            [-constraints.B; lpr], ...
                                    [], [], [], struct('info', 0));

end