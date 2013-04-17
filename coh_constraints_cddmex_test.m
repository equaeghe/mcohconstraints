function time = coh_constraints_cddmex_test(statenum, gamblenum)

  do
    try
      K = randomK(statenum, gamblenum);

      start = cputime;
        coh_constraints_cddmex(K);
      time = cputime - start;

      mistrial = 0;
    catch
      mistrial = 1;
    end_try_catch
  until mistrial == 0

end