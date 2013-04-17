function [time_vertrays, ...
          time_cddmex] = coh_constraints_cmp_test(statenum, gamblenum)

  do
    K = randomK(statenum, gamblenum);
    try
      start = cputime;
        coh_constraints_vertrays(K);
      time_vertrays = cputime - start;

      start = cputime;
        coh_constraints_cddmex(K);
      time_cddmex = cputime - start;

      mistrial = 0;
    catch
      mistrial = 1;
    end_try_catch
  until mistrial == 0

end