function coh_free_constraints_file(K, filename, numbertype)
% COH_FREE_CONSTRAINTS_FILE  writes the 'free' constraints for coherence to a file
%
% Synopsis:
%    coh_free_constraints_file(K, filename, numbertype)
%
% Input:
%    K = a nonnegative matrix with nonconstant columns ("gambles")
%    filename = the name of a file in the present working directory that will be
%               created (or overwritten!); suggested extension: .ext
%    numbertype = one of 'real', 'rational', or 'integer'; the number format
%                 used when writing the data to the file; the data is converted
%                 in case 'rational' or 'integer' is chosen, which in the latter
%                 case means 'truncated'
%
% Output:
%    a file named 'filename' in the present working directory
%
% Background & Method:
%    To be coherent, a lower prevision defined on a set of gambles K must
%    belong to the polyhedron defined by the constraints SP <= SK'μ_S,
%    μ_S >= 0, and 1'μ_S == 1 for all matrices S that differ from the identity
%    matrix by at most one signchange. These are written out to a file in
%    'Polyhedra H format'. This format can be read by polytpe theory programs
%    such as cddlib and lrs, for example, to compute the P-only
%    H-representation of the set of lower previsions on K that avoid sure loss.
%    Using cddlib, one would, e.g., say (on a Unix command line)
%
%           tail -1 my_output_file.ine | projection_gmp my_output_file.ine > my_projected_file.ine
%
%    (the last line of the output of this function is especially tailored
%    for this).

  [n, m] = size(K);
  N = 1:n;
  M = 1:m;

  % preallocate
  A = zeros((m + 1) * (m + n + 2), m * (1 + n));
  b = zeros((m + 1) * (m + n + 2), 1);


  A(M,M) = eye(m);
  A(M,m+N) = -K';
  A((m+1)*m+N,m+N) = -eye(n);
  A((m+1)*(m+n)+1,m+N) = ones(1, n);
  b((m+1)*(m+n)+1) = 1;
  A((m+1)*(m+n+1)+1,m+N) = -ones(1, n);
  b((m+1)*(m+n+1)+1) = -1;

  for k = M
    S = eye(m);
    S(k,k) = -1;

    A(k*m+M,M) = S;
    A(k*m+M,m+k*n+N) = -S*K';
    A((m+1)*m+k*n+N,m+k*n+N) = -eye(n);
    A((m+1)*(m+n)+k+1,m+k*n+N) = ones(1, n);
    b((m+1)*(m+n)+k+1) = 1;
    A((m+1)*(m+n+1)+k+1,m+k*n+N) = -ones(1, n);
    b((m+1)*(m+n+1)+k+1) = -1;
  end

  % depending on the output number type we need to convert the data and define
  % the format for writing out the data differently
  if strcmpi(numbertype, 'real')
    formatspec = [repmat(' % f', 1, m + (m + 1) * n + 1), '\n'];
  elseif strcmpi(numbertype, 'rational')
    A = rats(A);
    b = rats(b);
    formatspec = [repmat(' % s', 1, m + (m + 1) * n + 1), '\n'];
  elseif strcmpi(numbertype, 'integer')
    A = int32(A);
    b = int32(b);
    formatspec = [repmat(' % i', 1, m + (m + 1) * n + 1), '\n'];
  else
    error('invalid numbertype');
  end

  % open, write, and close the file
  fid = fopen(filename, 'wt');
    fprintf(fid, '%s\n', 'H-representation');
    fprintf(fid, '%s\n', 'begin');
    fprintf(fid, '%u %u %s\n', (m + 1) * (m + n + 2), ...
                               m + (m + 1) * n + 1, numbertype);
    fprintf(fid, formatspec, [b, -A]');
    fprintf(fid, '%s\n', 'end');
    fprintf(fid, ['%u ', repmat(' % u', 1, (m + 1) * n), '\n'], ...
                 (m + 1) * n, m+1:m+(m+1)*n);
  fclose(fid);

end