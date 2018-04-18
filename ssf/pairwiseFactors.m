function factors = pairwiseFactors(idxPairs,ES1,lambda1,ES2,lambda2,ES3,lambda3)
% converts all pairwise energies to factors

  nPairs   = size(idxPairs,1);
  factors  = cell(1,nPairs);

  for i=1:nPairs
      factors{i}.v = [idxPairs(i,1), idxPairs(i,2)];
      factors{i}.e = ...
        lambda1 * ES1{idxPairs(i,1), idxPairs(i,2)}(:)' + ...
        lambda2 * ES2{idxPairs(i,1), idxPairs(i,2)}(:)' + ...
        lambda3 * ES3{idxPairs(i,1), idxPairs(i,2)}(:)';

      if isempty(factors{i}.e)
          error('missing pairwise potentials')
      end
  end

end