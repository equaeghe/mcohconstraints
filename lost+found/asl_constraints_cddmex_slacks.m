function [lambdas, alphas] = asl_constraints_cddmex(K)
	% The asl-criterion (K-P)*lambda \not< 0 for all lambda >= 0
	% can be reduced to P*lambda <= alpha for all (lambda,alpha)-pairs
	% such that K*lambda + slacks = alpha such that lambda is pointwise
	% maximal and alpha is either 0 or 1. Slack variables are needed if
	% K is not of full rank, but as we don't test this here,
	% we always include them. For alpha=0, we need to bound lambda.

	[n, m] = size(K);

	% augment K to make sure we have a set of gambles of full rank
	KI = [K eye(n)]; % [! problematisch !]

	% positivity constraints x >= 0, i.e., -x <= 0
	x_geq_0 = [-eye(m+n) zeros(m+n, 1)];
	% boundedness constraints x <= 1
	x_leq_1 = [eye(m+n) ones(m+n, 1)];

	% case alpha=1: x >= 0 and KI*x = 1
	KI_x_eq_1 = [KI ones(n, 1)];
	% find vertices of feasible space A*x <= b
	Ab = [KI_x_eq_1; x_geq_0];
		% clear unneeded variables
		clear KI_x_eq_1;
	v1 = cddmex('extreme', struct('A', Ab(:,1:end-1), 'B', Ab(:,end), 'lin', (1:n)')).V;
	
	% case KI*x = 0 and alpha=0: 0 <= x <= 1
	KI_x_eq_0 = [KI zeros(n, 1)];
		% clear unneeded variables
		clear KI;
	% find vertices of feasible space A*x <= b
	Ab = [KI_x_eq_0; x_geq_0; x_leq_1];
		% clear unneeded variables
		clear x_geq_0 x_leq_1 KI_x_eq_0;
	v0 = cddmex('extreme', struct('A', Ab(:,1:end-1), 'B', Ab(:,end), 'lin', (1:n)')).V;

	% remove coefficients for the gambles that were added to guarantee full rank
	Ab = unique([v1(:,1:m) ones(size(v1, 1), 1); v0(:,1:m) zeros(size(v0, 1), 1)], 'rows');
		% clear unneeded variables
		clear v1 v0;

	% find miminal set of constraints [? add binary masks ?]
	h = cddmex('reduce_h', struct('A', Ab(:,1:end-1), 'B', Ab(:,end)));
	lambdas = h.A';
	alphas = h.B';
end