function K = randomK(statenum, gamblenum)

  Kevents = floor(rand(statenum, gamblenum) + rand);
  Kvals = floor(10 * rand(statenum, gamblenum));
  Kpos = Kevents .* Kvals;
  K = Kpos - ones(statenum, 1) * min(Kpos);

end