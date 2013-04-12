function constraints = coh_constraints_vertrays(K)
% COH_CONSTRAINTS_VERTRAYS  returns coherence constraints
%
% Synopsis:
%    coh_constraints_vertrays(K)
%
% Input:
%    K = a nonnegative matrix with nonconstant columns ("gambles")
%
% Output:
%    constraints = a struct describing constraints P'λ <= α with three
%                  fields:
%                    * A, a matrix containing the constraint coefficients λ
%                         as rows
%                    * B, a column vector containing the corresponding
%                         constraint constants α as components
%                    * lin, a column vector with indices of constraints
%                           that are actually equalities (== instead of <=)
%
% Background & Method:
%    To be coherent, a lower prevision defined on a set of gambles K must
%    for all matrices S that differ from the identity matrix by at most one
%    signchange be S-dominated by a linear prevision, i.e., a convex
%    combination of degenerate previsions; these correspond to the rows
%    of K. Here, S . This means that it belongs to the polyhedron defined
%    by the degenerate previsions as vertices and the columns of -S as
%    extreme rays.
%
%  See also CDDMEX, ASL_CONSTRAINTS_VERTRAYS

  [n, m] = size(K);

  % preallocate
  AB = cell(1,m);

  for k = 1:m
    % create S-matrix:
    S = eye(m);
    S(k,k) = -1;

    HS = cddmex('hull', struct('V', K, 'R', -S));

    ABS{k} = [HS.A, HS.B; -HS.A(HS.lin), -HS.B(HS.lin)];
  end

  AB = unique(vertcat(ABS{1:end}), 'rows');

  constraints = cddmex('reduce_h', ...
                       struct('A', AB(:,1:end-1), 'B', AB(:,end)));

end