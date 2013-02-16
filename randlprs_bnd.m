function lprs = randlprs_bnd(k, K)

  m = size(K, 2);

  Kmaxs = max(K)';
  Kmins = min(K)';
  Krngs = Kmaxs - Kmins;

  lprs = rand(m, k) .* repmat(Krngs, 1, k) - repmat(Kmins, 1, k);

end