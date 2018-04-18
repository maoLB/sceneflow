function sparseCrossTerm = crossTerm_flowDisp(S,D,F)
% project points to Ir{2} via left sparse flow observations + SGM 

  nS = numel(S);
  %[h,w] = size(D);

  Fu = F(:,:,1);
  Fv = F(:,:,2);
  Fi = F(:,:,3);

  sparseCrossTerm = [];

  for iS = 1:nS

    % project points to Il{2} via left flow
    
    % - find valid flow observations
    valid = Fi(S(iS).idx)>0;

    % - get flow from features
    fu = Fu(S(iS).idx(valid));
    fv = Fv(S(iS).idx(valid));

    % - project img coordinates
    uv2 = [S(iS).u(valid)+fu,S(iS).v(valid)+fv];
    uv2_int = round(uv2);

    u1 = S(iS).u(valid);
    v1 = S(iS).v(valid);

    % find valid disparity observations
    idx_d = sub2ind(size(D),uv2_int(:,2),uv2_int(:,1));
    valid = D(idx_d)>0;

    % avoid strange behaviour in case of only one invalid observation
    if sum(valid) == 0, continue; end

    % project points to Ir{2} via D{2}
    u4 = uv2(valid,1)-D(idx_d(valid));
    v4 = uv2(valid,2);

    sparseCrossTerm = [ sparseCrossTerm,  [u1(valid,:)';v1(valid,:)';u4';v4']  ];

  end

end