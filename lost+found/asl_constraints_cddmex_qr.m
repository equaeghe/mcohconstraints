function lambdas = asl_constraints_cddmex_qr(K)
	[L, Q1, Q2] = lq(K);
	mu = L \ ones(size(L, 1), 1);
	kappa = Q1' * mu;
	nus = cddmex('extreme', struct('A', -Q2', 'B', kappa));
	m = size(nus.V, 1);
	lambdas = repmat(kappa, 1, m) + Q2' * nus.V';
end