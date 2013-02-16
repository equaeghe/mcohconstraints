function lprs = randlprs_asl(k, K)

  [n, m] = size(K);

  mus = rand(n, k);
  munorms = sum(mus, 1);
  mus = mus ./ repmat(munorms, n, 1);

  lprs = K'*mus - rand(m, k);

end