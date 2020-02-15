function [basis, basis_indices, nonbasis_indices, ...
          nullspace] = matrixdigest_rref(A)
% matrixdigest_rref  returns a number of useful properties related matrices
%                    for a given matrix.
%
% Synopsis:
%    [basis, basis_indices, nonbasis_indices, ...
%     nullspace] = matrixdigest_rref(A)
%
% Input:
%    A = a matrix of any size containing numbers
%
% Output:
%    basis = a submatrix of A whose columns form a basis
%            for the column space of A
%    basis_indices = a row vector of indices to the columns of A that
%                    constitute the basis
%    nonbasis_indices = a row vector of indices to the columns of A not
%                       part of the basis
%    nullspace = a basis for the the (column/right) nullspace of A
%
% Remark:
%    This function makes use of Gauss-Jordan elimination under the hood.
%    For some matrices, this may introduce more round-off errors than other
%    methods would.

% Method:
%    We first bring A into reduced row echelon form. From this, we derive
%    the indices to a set of columns of A that form a basis for its column
%    space. The basis for the nullspace is then constructed from the
%    nonbasis columns of A.

  [n, m] = size(A);

  % To bring A into reduced row echelon form, we use the inbuilt function
  % 'rref'
  [reduced_A, basis_indices] = rref(A);

  rank = length(basis_indices);

  basis = A(:,basis_indices);

  nullity = m - rank;
  nonbasis_indices = 1:m;
    nonbasis_indices(basis_indices) = [];
  nullspace = zeros(m, nullity);
    nullspace(nonbasis_indices,:) = eye(nullity);
    nullspace(basis_indices,:) = -reduced_A(1:rank, nonbasis_indices);

end