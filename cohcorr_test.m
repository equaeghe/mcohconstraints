function cohcorr_test(fid, samplenum, statenum, gamblenum, lprnum)
% cohcorr_test  performs tests of the nadir point calculations
%
% Synopsis:
%    [Mtime0, Mtime, ...
%     Ftime0, Ftime] = cohcorr_test(statenum, gamblenum, lprnum)
%
% Input:
%    fid = file id, file to write output to 
%    statenum = an integer, the size of the possibility space
%    gamblenum = an integer, the size of the gamble space
%    lprnum = an integer, the numer of lower previsions for which a nadir
%             point must be calculated
%
% Output in file:
%    Mtime0 = the cputime needed to compute the minimal set of linear
%             constraints characterizing coherence
%    Mtime = the running time of the nadir calculations, divided by the
%            number of extreme points computed, averaged over the lprnum
%            trials, when using the minimal set of constraints
%    Ftime0 = the cputime needed to generate the 'free' set of linear
%             constraints characterizing coherence
%    Ftime = the running time of the nadir calculations, divided by the
%            number of extreme points computed, averaged over the lprnum
%            trials, when using the 'free' constraints
%    maximals = the number of maximal dominated coherent lower previsions
%
% Background & Method:
%
% See also COHCORR_BENSOLVE, RANDLPRS_BND,
%          COH_CONSTRAINTS_CDDMEX, COH_FREE_CONSTRAINTS.

  % how finely do we distinguish zero?
  eps = 10^(-8);

  k = 0;
  do

    % randomly select, without replacement, gamblenum NNZM gambles with
    % approximate sparsity .5
    K = randomK(statenum, gamblenum, .5);

    % calculate the minimal set of linear constraints characterizing
    % coherence
    start = cputime;
      Mconstraints = coh_constraints_vertrays(K);
    stop = cputime;
    Mtime0 = stop - start;

    % generate the 'free' set of linear constraints characterizing
    % coherence
    start = cputime;
      Fconstraints = coh_free_constraints(K);
    stop = cputime;
    Ftime0 = stop - start;

    % randomly generate a lower previsions
    lprs = randlprs_bnd(lprnum, K);

    % initial values
    Mtime = 0;
    Ftime = 0;

    % test the time it takes for a nadir computation for the lower
    % prevision generated
    for lpr = lprs

      try
        % when using the minimal set of constraints characterizing coherence
        start = cputime;
          [Mcohcorr, Mmaxnum] = cohcorr_bensolve(Mconstraints, lpr);
        stop = cputime;
        Mtime = Mtime + (stop - start);

        % when using the 'free' set of constraints characterizing coherence
        start = cputime;
          [Fcohcorr, Fmaxnum] = cohcorr_bensolve(Fconstraints, lpr);
        stop = cputime;
        Ftime = Ftime + (stop - start);

        excess = sum(abs(Mcohcorr - Fcohcorr));
        correction = sum(abs(Mcohcorr - lpr));
        if excess > eps
          % alert the user if both approaches give different nadir points
          K, lpr, Mcohcorr, Fcohcorr, excess
  %        error(['The minimal and free constraints give different' ...
  %               'answers (' num2str(Mmaxnum) ' vs. ' num2str(Fmaxnum) ')' ...
  %               'or roundoff error excess (' num2str(excess) ').']);
        elseif correction > eps
          % only record if the original lpr was incoherent
          fprintf(fid, '%1.1e\t%1.1e\t%1.1e\t%1.1e\t%u\r\n', ...
                  Mtime0, Mtime, Ftime0, Ftime, Mmaxnum);
          k += 1; k
        else
          correction
        end
      catch
        % drop cdd numerical instability-infected cases
      end_try_catch

    end

  until k >= lprnum * samplenum

end