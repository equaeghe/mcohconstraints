function [Mtime0, Mtime, ...
          Ftime0, Ftime] = cohcorr_test(statenum, gamblenum, lprnum)
% ***  ***
%
% Synopsis:
%    ***
%
% Input:
%    *** = ***
%
% Output:
%    *** = ***
%
% Background & Method:

  eps = 10^(-10);

  L = (dec2bin(1:2^statenum-2, statenum) - '0')';
  K = zeros(statenum, gamblenum);
  for k = 1:gamblenum
    l = ceil(rand * size(L, 2));
    K(:,k) = L(:,l);
    L(:,l) = [];
  end

  start = cputime;
    Mconstraints = coh_constraints_cddmex(K);
  stop = cputime;
  Mtime0 = stop - start;

  start = cputime;
    Fconstraints = coh_free_constraints(K);
  stop = cputime;
  Ftime0 = stop - start;

  lprs = randlprs_bnd(lprnum, K);

  Mtime = 0;
  Ftime = 0;

  for lpr = lprs

    start = cputime;
      [Mcohcorr, Mmaxnum] = cohcorr_bensolve(Mconstraints, lpr);
    stop = cputime;
    Mtime = Mtime + (stop - start) / Mmaxnum;

    start = cputime;
      [Fcohcorr, Fmaxnum] = cohcorr_bensolve(Fconstraints, lpr);
    stop = cputime;
    Ftime = Ftime + (stop - start) / Fmaxnum;

    if sum(abs(Mcohcorr - Fcohcorr)) > eps
      error(['The minimal and free constraints give different answers,' ...
             ' or roundoff error excess']);
    end
  end

  Mtime = Mtime / lprnum;
  Ftime = Ftime / lprnum;

end