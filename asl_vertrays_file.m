function asl_vertrays_file(K, filename, numbertype)
% ASL_VERTRAYS_FILE  writes out vertices and extreme rays for avoiding sure loss
%
% Synopsis:
%    asl_vertrays_file(K, filename, numbertype)
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
%    To avoid sure loss, a lower prevision defined on a set of gambles K
%    must be dominated by a linear prevision, i.e., a convex combination of
%    degenerate previsions; these correspond to the rows of K. This means
%    that it belongs to the polyhedron defined by the degenerate previsions
%    as vertices and negative directions as extreme rays. These are written
%    out to a file in 'Polyhedra V format'. This format can be read by
%    polytpe theory programs such as cddlib and lrs, for example, to
%    compute the facet or H-representation of the set of lower previsions
%    on K that avoid sure loss. Using cddlib, one would, e.g., say (on a
%    Unix command line)
%
%           lcdd_gmp my_asl_vertrays_file.ext > my_asl_vertrays_file.ine
%
%  See also ASL_CONSTRAINTS_VERTRAYS

  [n, m] = size(K);

  % put the data in one matrix, rows with leading ones indicate vertices,
  % rows with leading zeros indicate extreme rays
  data = [ones(n, 1) K; zeros(m, 1) -eye(m)];

  % depending on the output number type we need to convert the data and
  % define the format for writing out the data differently
  if strcmpi(numbertype, 'real')
    formatspec = [repmat(' % f', 1, m + 1), '\n'];
  elseif strcmpi(numbertype, 'rational')
    data = rats(data);
    formatspec = [repmat(' % s', 1, m + 1), '\n'];
  elseif strcmpi(numbertype, 'integer')
    data = int32(data);
    formatspec = [repmat(' % i', 1, m + 1), '\n'];
  else
    error('invalid numbertype');
  end

  % open, write, and close the file
  fid = fopen(filename, 'wt');
    fprintf(fid, '%s\n', 'V-representation');
    fprintf(fid, '%s\n', 'begin');
    fprintf(fid, '%u %u %s\n', n + m, m + 1, numbertype);
    fprintf(fid, formatspec, data');
    fprintf(fid, '%s\n', 'end');
  fclose(fid);

end