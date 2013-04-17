function time = coh_constraints_chain_test(statenum, gamblenum)

  do
    K = randomK(statenum, gamblenum);
    try
      start = cputime;
        coh_constraints_vertrays(K);
      time = cputime - start;
      mistrial = 0;
    catch
      try
        start = cputime;
          coh_constraints_cddmex(K);
        time = cputime - start;
        mistrial = 0;
      catch
        mistrial = 1;
      end_try_catch
    end_try_catch
  until mistrial == 0

end