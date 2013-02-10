function maximals = get_maximal_rows(unknowns)
% GET_MAXIMAL_ROWS  returns a matrix with pointwise maximal rows
%
% Synopsis:
%    maximals = get_maximal_rows(unknowns)
%
% Input:
%    unknowns = a matrix with rows of unknown relative pointwise
%               dominance status
%
% Output:
%    maximals = a matrix with the unique pointwise maximal rows of the
%               input matrix unknowns; relative order is not preserved
%
% Background & Method:
%    A vector x pointwise dominates another vector y if every component of
%    their difference x-y is nonpositive. This function uses a simple
%    scanning algorithm inspired by, but not equal to, the 'Best' algorithm
%    as presented in
%
%      Godfrey, P., Shipley, R., & Gryz, J. (2007). “Algorithms and analyses
%      for maximal vector computation.” The VLDB Journal, 16(1), 5-28.
%      doi:10.1007/s00778-006-0029-7

  % preallocate
  maximals = [];

  while ~isempty(unknowns)
    candidate = unknowns(1, :);

    % assume the candidate is maximal and compare it with all other rows
    maximals = [maximals; candidate];
    unknowns(1, :) = [];
    dominated = [];
    for k = 1:size(unknowns, 1)
      comparandum = unknowns(k,:);
      criterion = candidate - comparandum;
      % detect rows dominated by the candidate
      if all(criterion >= 0)
        dominated = [dominated; k];
      % detect if candidate is nonmaximal
      elseif all(criterion <= 0)
        maximals(end,:) = [];
        % "promote the winner"-heuristic: move dominating vector to the
        % front, to become the next candidate
        unknowns(k,:) = [];
        unknowns = [comparandum; unknowns];
        dominated = dominated + 1;
        break;
      end
    end
    % remove dominated from unknowns
    unknowns(dominated,:) = [];
  end

end