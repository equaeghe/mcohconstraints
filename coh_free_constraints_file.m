function coh_free_constraints_file(K, filename, numbertype)
% COH_FREE_CONSTRAINTS_FILE  writes out the 'free' coherence constraints
%
% Synopsis:
%    coh_free_constraints_file(K, filename, numbertype)
%
% Input:
%    K = a nonnegative matrix with nonconstant columns ("gambles")
%    filename = the name of a file in the present working directory that
%               will be created (or overwritten!);
%               suggested extension: .ext
%    numbertype = one of 'real', 'rational', or 'integer'; the number
%                 format used when writing the data to the file; the data
%                 is converted in case 'rational' or 'integer' is chosen,
%                 which in the latter case means 'truncated'
%
% Output:
%    a file named 'filename' in the present working directory
%
% Background & Method:
%    To be coherent, a lower prevision defined on a set of gambles K must
%    belong to the polyhedron defined by the constraints SP <= SK'μ_S,
%    μ_S >= 0, and 1'μ_S == 1 for all matrices S that differ from the
%    identity matrix by at most one signchange. These are written out to a
%    file in 'Polyhedra H format'. This format can be read by polytpe
%    theory programs such as cddlib and lrs, for example, to compute the
%    P-only H-representation of the set of lower previsions on K that avoid
%    sure loss. Using cddlib, one would, e.g., say (on a Unix command line)
%
%       tail -1 my_output_file.ine | projection_gmp my_output_file.ine \
%                                  > my_projected_file.ine
%
%    (the last line of the output of this function is especially tailored
%    for this).
%
% See also COH_FREE_CONSTRAINTS, ASL_FREE_CONSTRAINTS_FILE

  constraints = asl_free_constraints(K);
  A = constraints.A;
  b = constraints.B;
  clear constraints;
  [n, m] = size(K);
  [Arows, Acols] = size(A);

  % depending on the output number type we need to convert the data and
  % define the format for writing out the data differently
  if strcmpi(numbertype, 'real')
    formatspec = [repmat(' % f', 1, Acols + 1), '\n'];
  elseif strcmpi(numbertype, 'rational')
    A = rats(A);
    b = rats(b);
    formatspec = [repmat(' % s', 1, Acols + 1), '\n'];
  elseif strcmpi(numbertype, 'integer')
    A = int32(A);
    b = int32(b);
    formatspec = [repmat(' % i', 1, Acols + 1), '\n'];
  else
    error('invalid numbertype');
  end

  % open, write, and close the file
  fid = fopen(filename, 'wt');
    fprintf(fid, '%s\n', 'H-representation');
    fprintf(fid, '%s\n', 'begin');
    fprintf(fid, '%u %u %s\n', Arows, Acols + 1, numbertype);
    fprintf(fid, formatspec, [b, -A]');
    fprintf(fid, '%s\n', 'end');
    fprintf(fid, ['%u ', repmat(' % u', 1, (m + 1) * n), '\n'], ...
                 (m + 1) * n, m+1:m+(m+1)*n);
  fclose(fid);

end