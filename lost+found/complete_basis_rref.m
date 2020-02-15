function B = complete_basis_rref(A)
% complete_basis_rref  returns a matrix with column vectors that are
%                      orthogonal to the columns of the given matrix.
%
% Synopsis:
%    B = complete_basis_rref(A)
%
% Input:
%    A = a matrix of any size containing numbers
%
% Output:
%    B = a matrix whose columns are orthogonal to the columns of A and
%        such that [A B] has full rank
%
% Remark:
%    This function makes use of Gauss-Jordan elimination under the hood.
%    For some matrices, this may introduce more round-off errors than other
%    methods would.

% Method:
%    The row nullspace of a matrix A is orthogonal to its column space,
%    so to complete a basis, we can use a basis for this nullspace.
%    This basis is calculated here by first bringing A into reduced row
%    echelon form.

  [~, ~, ~, ~, ~, ~, B] = matrixdigest_rref(A');

end