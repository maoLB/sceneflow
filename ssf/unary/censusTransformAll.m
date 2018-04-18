function [ Cl,Cr,time_spent ] = censusTransformAll( Il,Ir )
%CENSUSTRANSFORMALL Computes census transformed images for four views

tStart = tic;

[~,~,ch] = size(Il{1});
assert(ch==1,'input images must be single-channel');

n = numel(Il);

Cl = cell(1,n);
Cr = cell(1,n);

for frame=1:n
   Cl{frame} = colfilt(Il{frame},[5 5],'sliding',@census);
   Cr{frame} = colfilt(Ir{frame},[5 5],'sliding',@census);
end

time_spent = toc(tStart);
% fprintf('\nCensus transform computed in\t%6.2f s\n',time_spent);

end

