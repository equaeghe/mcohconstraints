function K = randomK(statenum, gamblenum, cutoff)

  Kevents = floor(rand(statenum, gamblenum) + cutoff);
  Kvals = floor(10 * rand(statenum, gamblenum));
  Kpos = Kevents .* Kvals;
  K = Kpos - ones(statenum, 1) * min(Kpos);

end