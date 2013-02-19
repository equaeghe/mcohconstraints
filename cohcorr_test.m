function [Mtime0, Mtime, ...
          Ftime0, Ftime] = cohcorr_test(statenum, gamblenum, lprnum)
% cohcorr_test  performs tests of the nadir point calculations
%
% Synopsis:
%    [Mtime0, Mtime, ...
%     Ftime0, Ftime] = cohcorr_test(statenum, gamblenum, lprnum)
%
% Input:
%    statenum = an integer, the size of the possibility space
%    gamblenum = an integer, the size of the gamble space
%    lprnum = an integer, the numer of lower previsions for which a nadir
%             point must be calculated
%
% Output:
%    Mtime0 = the cputime needed to compute the minimal set of linear
%             constraints characterizing coherence
%    Mtime = the running time of the nadir calculations, divided by the
%            number of extreme points computed, averaged over the lprnum
%            trials, when using the minimal set of constraints
%    Ftime0 = the cputime needed to generate the 'free' set of linear
%             constraints characterizing coherence
%    Mtime = the running time of the nadir calculations, divided by the
%            number of extreme points computed, averaged over the lprnum
%            trials, when using the 'free' constraints
%
% Background & Method:
%
% See also COHCORR_BENSOLVE, RANDLPRS_BND,
%          COH_CONSTRAINTS_CDDMEX, COH_FREE_CONSTRAINTS.

  % how finely do we distinguish zero?
  eps = 10^(-10);

  % randomly select, without replacement, gamblenum indicator functions
  L = (dec2bin(1:2^statenum-2, statenum) - '0')';
  K = zeros(statenum, gamblenum);
  for k = 1:gamblenum
    l = ceil(rand * size(L, 2));
    K(:,k) = L(:,l);
    L(:,l) = [];
  end

  % calculate the minimal set of linear constraints characterizing
  % coherence
  start = cputime;
    Mconstraints = coh_constraints_cddmex(K);
  stop = cputime;
  Mtime0 = stop - start;

  % generate the 'free' set of linear constraints characterizing
  % coherence
  start = cputime;
    Fconstraints = coh_free_constraints(K);
  stop = cputime;
  Ftime0 = stop - start;

  % randomly generate the lower previsions
  lprs = randlprs_bnd(lprnum, K);

  % initial values
  Mtime = 0;
  Ftime = 0;

  % test the time it takes for a nadir computation for each lower
  % previsions generated
  for lpr = lprs

    % when using the minimal set of constraints characterizing coherence
    start = cputime;
      [Mcohcorr, Mmaxnum] = cohcorr_bensolve(Mconstraints, lpr);
    stop = cputime;
    Mtime = Mtime + (stop - start) / Mmaxnum;

    % when using the 'free' set of constraints characterizing coherence
    start = cputime;
      [Fcohcorr, Fmaxnum] = cohcorr_bensolve(Fconstraints, lpr);
    stop = cputime;
    Ftime = Ftime + (stop - start) / Fmaxnum;

    % alert the user if both approaches give different nadir points
    if sum(abs(Mcohcorr - Fcohcorr)) > eps
      error(['The minimal and free constraints give different answers,' ...
             ' or roundoff error excess']);
    end
  end

  % compute average
  Mtime = Mtime / lprnum;
  Ftime = Ftime / lprnum;

end