function [L, Q1, Q2] = lq(A)
	n = size(A, 1);
	[Q, R] = qr(A');
	L = R(1:n,:)';
	Q1 = Q(:, 1:n)';
	Q2 = Q(:, n+1:end)';
end