function [lambdas, alphas] = coh_constraints_cddmex(K)
	[n, m] = size(K);

	% augment K to make sure we have a set of gambles of full rank
	KI = [K eye(n)];

	% positivity constraints x >= 0, i.e., -x <= 0
	x_geq_0 = [-eye(m+n) zeros(m+n, 1)];
	% boundedness constraints x <= 1
	x_leq_1 = [eye(m+n) ones(m+n, 1)];

	for k = 1:m

		% case KIk*x = 1
		% positivity constraints x >= 0, i.e., -x <= 0
		x_geq_0k = [-eye(m+n) zeros(m+n, 1)];
		KIk = KI;
		KIk(:,k) = -KIk(:,k);
		KIk_x_eq_1 = [KIk ones(n, 1)];
			% clear unneeded variables
			clear KIk;
		% find vertices of feasible space A*x <= b
		Ab = [KIk_x_eq_1; x_geq_0k];
			% clear unneeded variables
			clear KI_x_eq_1 x_geq_0k;
		v1 = cddmex('extreme', struct('A', Ab(:,1:end-1), 'B', Ab(:,end), 'lin', (1:n)')).V;
		% invert f coefficient, remove rank-filling coefficients
		lambdas1{k} = unique([v1(:,1:k-1), -v1(:,k), v1(:,k+1:m)], 'rows')';
			% clear unneeded variables
			clear v1;

		% case KIk*x = f
		% positivity constraints x >= 0, i.e., -x <= 0
		x_geq_0k = [-eye(m-1+n) zeros(m-1+n, 1)];
		KIk = KI;
		KIk(:,k) = [];
		KIk_x_eq_f = [KIk K(:,k)];
			% clear unneeded variables
			clear KIk;
		% find vertices of feasible space A*x <= b
		Ab = [KIk_x_eq_f; x_geq_0k];
			% clear unneeded variables
			clear KI_x_eq_f x_geq_0k;
		v0 = cddmex('extreme', struct('A', Ab(:,1:end-1), 'B', Ab(:,end), 'lin', (1:n)')).V;
		% reinsert f coefficient, remove rank-filling coefficients
		lambdas0{k} = unique([v0(:,1:k-1), -ones(size(v0, 1), 1), v0(:,k:m-1)], 'rows')';
			% clear unneeded variables
			clear v0;

		% case KIk*x = -1
		% positivity constraints x >= 0, i.e., -x <= 0
		x_geq_0k = [-eye(m+n) zeros(m+n, 1)];
		KIk = KI;
		KIk(:,k) = -KIk(:,k);
		KIk_x_eq_min1 = [KIk -ones(n, 1)];
			% clear unneeded variables
			clear KIk;
		% find vertices of feasible space A*x <= b
		Ab = [KIk_x_eq_min1; x_geq_0k];
			% clear unneeded variables
			clear KI_x_eq_min1 x_geq_0k;
		vmin1 = cddmex('extreme', struct('A', Ab(:,1:end-1), 'B', Ab(:,end), 'lin', (1:n)')).V;
		% invert f coefficient, remove rank-filling coefficients
		lambdasmin1{k} = unique([vmin1(:,1:k-1), -vmin1(:,k), vmin1(:,k+1:m)], 'rows')';
			% clear unneeded variables
			clear vmin1;

	end

	lambdas1 = [lambdas1{1:end}]
	lambdas0 = [lambdas0{1:end}]
	lambdasmin1 = [lambdasmin1{1:end}]

	Ab = [lambdas1' ones(size(lambdas1, 2), 1); lambdas0' ones(size(lambdas0, 2), 1); lambdasmin1' -ones(size(lambdasmin1, 2), 1)];

	% positivity constraints lambda >= 0, i.e., -lambda <= 0
%	lambda_geq_0 = [-eye(m) zeros(m, 1)];

	% find miminal set of constraints, taking positivity into account
%	Ab = [Ab; lambda_geq_0];
	h = cddmex('reduce_h', struct('A', Ab(:,1:end-1), 'B', Ab(:,end)));
	coh_lambdas = h.A';
	coh_alphas = h.B';

	[asl_lambdas, asl_alphas] = asl_constraints_cddmex(K);

	h = cddmex('reduce_h', struct('A', [asl_lambdas coh_lambdas]', 'B', [asl_alphas coh_alphas]'));
	lambdas = h.A';
	alphas = h.B';

end
