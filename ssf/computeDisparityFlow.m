function [ dD,D2,S2 ] = computeDisparityFlow( S,shape,object,h,w,K,base_m )
%COMPUTEDISPARITYFLOW computes disparity flow map and warped disp map at t1

  S2 = S;
  L2 = zeros(h,w);
  nS = numel(S2);

  % init disp map
  D1      = zeros(h,w);
  D2      = zeros(h,w);
  D2_noc  = zeros(h,w);


  for iS = 1:nS % all superpixels

    % disparities at t0
    u1 = S(iS).u;
    v1 = S(iS).v;   

    % compute all disparities
    d  = shape(iS,1) * (u1-S(iS).cu) + shape(iS,2) * (v1-S(iS).cv) + shape(iS,3);

    % assign disparities
    D1(S(iS).idx) = d; 

    % get 3D normal at t0
    dispPlane_abc = abcFromAlphaBetaGamma(shape(iS,1:3),[S(iS).cu,S(iS).cv]);
    n{1} = dispPlaneToNormal(dispPlane_abc,K,base_m);

    % get transformation matrix from object motion
    k = shape(iS,4);

    rx_d    = object(k,1);
    ry_d    = object(k,2);
    rz_d    = object(k,3);
    tx      = object(k,4);
    ty      = object(k,5);
    tz      = object(k,6);

    Tr = getRigidMotionTrafo(rx_d,ry_d,rz_d,tx,ty,tz);

    R  = Tr(1:3,1:3);
    t  = Tr(1:3,4);

    % transform the plane to t2
    d1 = 1/norm(n{1});
    n1 = n{1}/norm(n{1});
    p1 = -d1.*n1;
    p2 = Tr*[p1;1];
    n2 = Tr*[n1;0];
    d2 = -(n2'*p2);

    n{2} = n2(1:3)./d2;

    abc = normalToDispPlane(n{2},K,base_m);

    % compute all disparities
    d2  = abc(1)*u1+abc(2)*v1+abc(3);

    % save disparity map at t1 induced by transformed planes
    D2(S(iS).idx) = d2; 

    % warp image coordinates
    H  = (K * (R-t*n{1}') / K);

    uvh = H * [u1';v1';ones(1,length(u1))];
    u2 = round((uvh(1,:)./uvh(3,:))');
    v2 = round((uvh(2,:)./uvh(3,:))');

    d2 = abc(1)*u2+abc(2)*v2+abc(3);

    % assign disparities
    D2_noc(S(iS).idx) = max(D2_noc(S(iS).idx),d2); 

    % remove points outside the image plane
    invalid = v2<1 | v2>h | u2<1 | u2>w;

    u2( invalid ) = [];
    v2( invalid ) = [];

    idx2 = sub2ind(size(D2_noc),v2,u2);

    % warp labelling
    L2(idx2) = iS;

    % set parameters of transformed superpixel
    S2(iS).n  = n{2};
    S2(iS).nP = numel(u2);
    S2(iS).cu = median(u2);
    S2(iS).cv = median(v2);

    alphaBetaGamma = alphaBetaGammaFromAbc(abc,[S2(iS).cu,S2(iS).cv]);
    S2(iS).alpha   = alphaBetaGamma(1);
    S2(iS).beta    = alphaBetaGamma(2);
    S2(iS).gamma   = alphaBetaGamma(3);

  end

  % compute disparity flow
  dD = D2_noc-D1;

end

