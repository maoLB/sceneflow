function [F,motionMap] = flowMapFromShapeAndMotion( S,shape,object,h,w,K,base_m )
%FLOWMAPFROMSHAPEANDMOTION Flow map from superpixels,
% (alpha,beta,gamma)-planes and motion particles

nS     = numel(S);

% init disp and label map
Fu = zeros(h,w);
Fv = zeros(h,w);

motionMap = zeros(h,w);

planes_abc = abcFromAlphaBetaGamma( shape(:,1:3), [S.cu;S.cv]' );

for iS = 1:nS  % all superpixels
    
    % index of optimal motion
    iM = shape(iS,4);
    motionMap(S(iS).idx) = iM;
    
    % get the 3D plane normal
    normal = dispPlaneToNormal( planes_abc(iS,:),K,base_m );

    % get transformation matrix from motion parameters
    rx_d    = object(iM,1,end);
    ry_d    = object(iM,2,end);
    rz_d    = object(iM,3,end);
    tx      = object(iM,4,end);
    ty      = object(iM,5,end);
    tz      = object(iM,6,end);

    Tr = getRigidMotionTrafo(rx_d,ry_d,rz_d,tx,ty,tz);

    % compute homography induced by assigned normal+motion
    Rc = Tr(1:3,1:3);
    tc = Tr(1:3,4);
    
    H  = (K * (Rc-tc*normal') / K);
          
    u1 = S(iS).u;
    v1 = S(iS).v;
    
    uvh = H * [u1';v1';ones(1,length(u1))];
    u2 = ((uvh(1,:)./uvh(3,:))');
    v2 = ((uvh(2,:)./uvh(3,:))');
    
    Fu(S(iS).idx) = u2-u1; 
    Fv(S(iS).idx) = v2-v1; 
end

F = cat(3,Fu,Fv,ones(size(Fu)));

end

