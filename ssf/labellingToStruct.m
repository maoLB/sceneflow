function [S,Dp,distS] = labellingToStruct(L,D)
% computes mean disparity & dispPlane parameters from L & D
% returns struct S containing parameterized superpixels

% number of superpixels
nS = max(L(:));
[h,w] = size(L);

% init superpixel struct and disparity map
s  = struct();
S(1:nS) = s;
Dp = zeros(h,w);
valid = true(nS,1);

% compute gradients of the disparity map
% pad result with zeros since diff reduces one dimension by 1
D_dv = [diff(D,1,1);zeros(1,w)];
D_du = [diff(D,1,2),zeros(h,1)];

% for all superpixels do
for iS=1:nS
  
    % save indices of member pixels
    idx = find(L==iS);
    S(iS).idx = idx;

    % get all pixels belonging to current superpixel
    [S(iS).v,S(iS).u] = ind2sub(size(L),idx);

    % set number of pixels
    S(iS).nPix  = length(idx);

    % set median u,v as approximation of superpixel center
    S(iS).cu = median(S(iS).u);
    S(iS).cv = median(S(iS).v);
    
    % get disparities from input disp map
    d = D(idx);
    
    % boolean indicator of valid disparities
    % assuming invalid disparities are marked as -1 / negative values
    iV = d >= 0;
   
    % if more than 50% of the pixels have a valid disparity get disparity
    % plane parameters
    if sum(iV) > 0.5*numel(S(iS).idx) 

        % alpha, beta, gamma refer to the centered representation of 
        % the disparity plane as used in Yamaguchi's 2012 ECCV paper
        S(iS).alpha   = median(D_du(idx));
        S(iS).beta    = median(D_dv(idx));
        S(iS).gamma   = median(d - (S(iS).alpha * (S(iS).u - S(iS).cu) + S(iS).beta * (S(iS).v - S(iS).cv)));
        Dp(idx) = S(iS).gamma + S(iS).alpha * (S(iS).u - S(iS).cu) + S(iS).beta * (S(iS).v - S(iS).cv);
        
    else % mark the disparity plane as invalid
        valid(iS)   = false;
        S(iS).alpha =  0;
        S(iS).beta  =  0;
        S(iS).gamma = -1;
    end

end

% init invalid superpixels with the disparity of its nearest neighbor
centers     = [S.cu;S.cv]';
gammas      = [S.gamma];
validGammas = gammas(valid);
nnIdx       = knnsearch(centers(valid,:), centers(~valid,:));

c = 1;
for iS = 1:nS
    if ~valid(iS)
        S(iS).gamma = validGammas(nnIdx(c));
        Dp(S(iS).idx) = S(iS).gamma + S(iS).alpha * S(iS).u + S(iS).beta * S(iS).v;
        c = c+1;
    end
end

% compute distances between all pairs of superpixels
distS = pdist2(centers,centers);
