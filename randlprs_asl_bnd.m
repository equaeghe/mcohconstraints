function lprs = randlprs_asl_bnd(k, K)

  Kmins = min(K)';

  lprs = randlprs_asl(k, K);
  lprs = max(lprs, repmat(Kmins, 1, k));

end