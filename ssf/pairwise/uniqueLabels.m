function [ L ] = uniqueLabels( L )

  % get unique labels
  uL = unique(L{1}(:));
  
  % get number of unique labels
  nL = double(numel(uL));
  
  l=zeros(size(L{1}));
  
  for i=1:nL
    l(L{1}==uL(i))=i;
  end
  
  L{1}=l;
  
end

