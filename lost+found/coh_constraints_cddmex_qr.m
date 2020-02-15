function [C1, C0] = coh_constraints_cddmex_qr(K)
	[cardinality, k] = size(K);

	% RHS=1-case
	C1 = [];
	[L, Q1, Q2] = lq(K);
	mu1 = L \ ones(cardinality, 1);
	for n = 1:k
		Q1n = Q1;
		Q1n(:,n) = -Q1(:,n);
		Q2n = Q2;
		Q2n(:,n) = -Q2(:,n);
		kappa = Q1n' * mu1;
		nus = cddmex('extreme', struct('A', -Q2n', 'B', kappa));
		m = size(nus.V, 1);
		Cn = repmat(kappa, 1, m) + Q2n' * nus.V';
		Cnn = [Cn(1:n-1,:); -Cn(n,:); Cn(n+1:k,:)];
		C1 = [C1, Cnn];
	end

	% RHS=0-case
	C0 = [];
	for n = 1:k
		Kn = K;
		Kn(:,n) = [];
		[Ln, Q1n, Q2n] = lq(Kn);
		mun = Ln \ K(:,n);
		kappan = Q1n' * mun;
		nus = cddmex('extreme', struct('A', -Q2n', 'B', kappan));
		m = size(nus.V, 1);
		Cn = repmat(kappan, 1, m) + Q2n' * nus.V';
		Cnn = [Cn(1:n-1,:); -ones(1,m); Cn(n:k-1,:)];
		C0 = [C0, Cnn];
	end
end
