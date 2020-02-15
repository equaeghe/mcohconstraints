function [lambdas, alphas] = coh_constraints_cddmex_rref(K)

  % calculate the number of rows and columns, the rank, a basis, indices to
  % the basis vectors, and the nullspace of K
  [k, m, r, B, b, N] = matrixdigest_rref(K);

  % if K does not have full rank, augment it so that it does
  if r < k
    [~, m, ~, B, b, N] = matrixdigest_rref([complete_basis_rref(B), K]);
  end

  for n = 1:m+1
    % 
    Bn = B;
    Nn = N;
    if n ~= m+1
      if any(n == b)
        Bn(:,n) = -B(:,n);
        Nn(n,:) = -N(n,:);
      else
        Nn(b,n - k) = -N(b,n - k);
      end
    end

    % 0-case
    kappa0 = zeros(m, 1);
    mus0 = cddmex('extreme', struct('A', [-N; N], ...
                                    'B', [kappa0; ones(m, 1)])).V;
    lambdas0{n} = mus0 * N';
    lambdas0{n}(:,n) = -lambdas0{n}(:,n);
    clear kappa0 mus0;

    % 1-case
    kappa1 = zeros(m, 1);
    kappa1(b) = B \ ones(k, 1);
    mus1 = cddmex('extreme', struct('A', -N, 'B', kappa1)).V;
    lambdas1{n} = repmat(kappa1', size(mus1, 1), 1) + mus1 * N';
    lambdas1{n}(:,n) = -lambdas1{n}(:,n);
    clear kappa1 mus1;

    % -1-case
    kappa2 = zeros(m, 1);
    kappa2(b) = B \ -ones(k, 1);
    mus2 = cddmex('extreme', struct('A', -N, 'B', kappa2)).V;
    lambdas2{n} = repmat(kappa2', size(mus2, 1), 1) + mus2 * N';
    lambdas2{n}(:,n) = -lambdas2{n}(:,n);
    clear kappa2 mus2;

    if r < k
      lambdas0{n} = unique(lambdas0{n}(:,k-r+1:end), 'rows');
      lambdas1{n} = unique(lambdas1{n}(:,k-r+1:end), 'rows');
      lambdas2{n} = unique(lambdas2{n}(:,k-r+1:end), 'rows');
    end
  end

  lambdas0 = unique(vertcat(lambdas0{1:end}), 'rows');
  lambdas1 = unique(vertcat(lambdas1{1:end}), 'rows');
  lambdas2 = unique(vertcat(lambdas2{1:end}), 'rows');

  h = cddmex('reduce_h', struct('A', [lambdas0; lambdas1;
                                      lambdas2; -eye(m)], ...
                                'B', [zeros(size(lambdas0, 1), 1); ...
                                      ones(size(lambdas1, 1), 1); ...
                                      -ones(size(lambdas2, 1), 1); ...
                                      zeros(m, 1)]));
  clear lambdas0 lambdas1 lambdas2;
  lambdas = h.A';
  alphas = h.B';
end