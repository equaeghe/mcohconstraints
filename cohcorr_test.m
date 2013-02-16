function [Mtime, Ftime] = cohcorr_test(statenum, gamblenum, lprnum)

  eps = 10^(-10);

  L = (dec2bin(1:2^statenum-2, statenum) - '0')';
  K = zeros(statenum, gamblenum);
  for k = 1:gamblenum
    l = ceil(rand * size(L, 2));
    K(:,k) = L(:,l);
    L(:,l) = [];
  end

  Mconstraints = coh_constraints_cddmex(K);
  Fconstraints = coh_free_constraints(K);

  lprs = randlprs_bnd(lprnum, K);

  Mtime = 0;
  Ftime = 0;

  for lpr = lprs

    Mtime = Mtime - cputime;
      Mcohcorr = cohcorr_bensolve(Mconstraints, lpr);
    Mtime = Mtime + cputime;

    Ftime = Ftime - cputime;
      Fcohcorr = cohcorr_bensolve(Mconstraints, lpr);
    Ftime = Ftime + cputime;

    if sum(abs(Mcohcorr - Fcohcorr)) > eps
      error(['The minimal and free constraints give different answers,' ...
             ' or roundoff error excess']);
    end
  end

end