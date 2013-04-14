function [time_vertrays, ...
          time_cddmex, ...
          mistrials] = coh_constraints_cmp_test(statenum, gamblenum, trials)

  % initialize
  time_vertrays = 0;
  time_cddmex = 0;
%  time_bensolve = 0;

  k = 0;
  mistrials = 0;
  do
    try
      % generate random K
      Kevents = floor(rand(statenum, gamblenum) + rand);
      Kvals = floor(10 * rand(statenum, gamblenum));
      Kpos = Kevents .* Kvals;
      K = Kpos - ones(statenum, 1) * min(Kpos);

      start = cputime;
        coh_constraints_vertrays(K);
      stop = cputime;
      time_vertrays = time_vertrays + stop - start;

      start = cputime;
        coh_constraints_cddmex(K);
      stop = cputime;
      time_cddmex = time_cddmex + stop - start;

  %    start = cputime;
  %      coh_constraints_bensolve(K);
  %    stop = cputime;
  %    time_bensolve = time_bensolve + stop - start;
      k += 1;
    catch
      mistrials += 1;
    end_try_catch
  until k == trials

  time_vertrays = time_vertrays / trials;
  time_cddmex = time_cddmex / trials;
%  time_bensolve = time_bensolve / trials;

end