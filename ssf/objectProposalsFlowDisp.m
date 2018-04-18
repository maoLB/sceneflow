function [Tr,objects,time_spent] = ...
  objectProposalsFlowDisp(sparseFlowLeft,D,~,F,Pd,s)
%objectProposalsFlowDisp generates object hypotheses
% generates hypotheses of individually moving objects from consistently
% moving clusters of sparse flow features that don't conform to bg motion

  t1 = tic;

	%% Parameters
  dispThreshLo    =   5;
  flowThresh      =   2;
  radiusCluster   =   2.5; % 3D radius of points considered for Ransac
  threshNonmax    =   3;
  
  nHyp            =  50;   % number of object hypothesis seed points
    
  % Ransac parameters
  nRansac  = 200; % number of RANSAC iterations
  threshIn = 2;   % inlier threshold RANSAC [px]
  nIn      = 5;   % number of inliers required to generate a proposal
    
  %% Get valid disparities
  % get disparities
  uv1  = sparseFlowLeft(1:2,:)';
  uv2  = sparseFlowLeft(3:4,:)';
  idx1 = sub2ind(size(D{1}),uv1(:,2),uv1(:,1));
  idx2 = sub2ind(size(D{1}),uv2(:,2),uv2(:,1));

  d1   = D{1}(idx1);
  d2   = D{2}(idx2);

  % find points with valid disparites
  validDisp = d1 > dispThreshLo & d2 > dispThreshLo;

  uv1  = uv1(validDisp,:);
  uv2  = uv2(validDisp,:);

  idx1 = idx1(validDisp);

  d1   = d1(validDisp); 
  d2   = d2(validDisp);
    
  %% Remove flow vectors that are explained by background motion
  fu = F(:,:,1);    
  fv = F(:,:,2);

  % compute flow components
  fu_sp = uv2(:,1)-uv1(:,1);
  fv_sp = uv2(:,2)-uv1(:,2);

  % compute optical flow endpoint error
  du = fu_sp-fu(idx1);
  dv = fv_sp-fv(idx1);
  epe = sqrt( (du).^2 + (dv).^2 );

  % find points contradicting background motion
  validFlow = epe > flowThresh;

  uv1  = uv1(validFlow,:);
  uv2  = uv2(validFlow,:);

  d1   = d1(validFlow); 
  d2   = d2(validFlow);

  % project points to 3D
  X1_all = project([uv1,d1],inv(Pd));
  X2_all = project([uv2,d2],inv(Pd));

  nPoints = size(X1_all,1);
  
  Tr_ = [];
  i   = 1;
  nIter = 0;

  %% Handle all background case -> no outliers
  if nPoints < 10
    Tr = [];
    objects = [];
    time_spent = toc(t1);
    return;
  end

  %% RANSAC loop
  nInliers = zeros(nHyp,1);
  
  while nIter < nHyp
              
    nIter = nIter+1;
        
    % sample cluster center
    if exist('s','var') && ~isempty(s)
      idxCC = randi(s,nPoints);
    else
      idxCC = randi(nPoints);
    end
        
    % get points around center
    dCC   = sqrt( (X1_all(:,1)-X1_all(idxCC,1)).^2 ...
                + (X1_all(:,2)-X1_all(idxCC,2)).^2 ...
                + (X1_all(:,3)-X1_all(idxCC,3)).^2 );
                
    idxCC = find(dCC < radiusCluster);
        
    X1 = X1_all(idxCC,:);
    X2 = X2_all(idxCC,:);
        
    uv1_c = uv1(idxCC,:);
    uv2_c = uv2(idxCC,:);
        
    nPointsCC = size(X1,1);
        
    % continue if there are not enough points in the cluster
    if nPointsCC<3, continue; end

    nInliers(i) = nIn;
    success = 0;

    for j=1:nRansac

      % draw indices
      if exist('s','var') && ~isempty(s)
        idxRansac = randperm(s,nPointsCC,3);
      else
        idxRansac = randperm(nPointsCC);
        idxRansac = idxRansac(1:3);
      end

      % estimate rigid motion
      Tr_c = getRigidMotion(X1(idxRansac,:),X2(idxRansac,:));

      % find inliers
      X2_c  = project(X1,Tr_c);
      uvd_l = project(X2_c,Pd);

      % reprojection error in all 4 images (covers depth component)
      uvd2 = project(X2,Pd);
      du   = uvd2(:,1) - uvd_l(:,1);
      dv   = uvd2(:,2) - uvd_l(:,2);
      dd   = uvd2(:,3) - uvd_l(:,3);
      dil  = sqrt(du.^2+dv.^2+dd.^2);

      % 3D error
      dX   = X2(:,1)-X2_c(:,1);
      dY   = X2(:,2)-X2_c(:,2);
      dZ   = X2(:,3)-X2_c(:,3);
      d3d  = sqrt(dX.^2+dY.^2+dZ.^2);

      nI   = sum( dil<threshIn );

      if nI > nInliers(i) % valid hypothesis      
        
        nInliers(i) = nI;
        idxIn   = find( dil<threshIn );
        Tr_{i}  = Tr_c;
                
        % cluster characteristics
        uv1_cluster{i} = uv1_c(idxIn,:);
        uv2_cluster{i} = uv2_c(idxIn,:);
        X1_cluster{i}  = X1(idxIn,:);
        X2_cluster{i}  = X2(idxIn,:);
        X1_median(i,:) = median(X1(idxIn,:));
        error2D{i} = dil(idxIn);
        error3D{i} = d3d(idxIn);
        success = 1;
      end
    end
        
    if success==1, i=i+1; end
  end
    
  time_spent = toc(t1);
  %fprintf('Finished RANSAC in\t\t%6.2f s\n',time_spent);

  % if there is not a single object hypothesis
  % return empty matrices
  if isempty(Tr_)
    Tr = [];
    objects = [];
    return;    
  end
  
  %% Postprocessing of clusters
  t2 = tic;

  % sort clusters wrt. inliers
  [~,idxMaxIn] = sort(nInliers,'descend');
    
  % always accept the biggest cluster
  Tr{1}   = Tr_{idxMaxIn(1)};
  XZ(1,:) = [ X1_median(idxMaxIn(1),1), X1_median(idxMaxIn(1),3) ];
  uv1_max{1} = uv1_cluster{idxMaxIn(1)};
  uv2_max{1} = uv2_cluster{idxMaxIn(1)};

  objects{1}.Tr = Tr_{idxMaxIn(1)};
  objects{1}.uv = uv1_cluster{idxMaxIn(1)};
  objects{1}.X1 = X1_cluster{idxMaxIn(1)};
  objects{1}.X2 = X2_cluster{idxMaxIn(1)};
  objects{1}.X_median = X1_median(idxMaxIn(1),:);
    
  % add further clusters using nonmax-suppression
  iTr = 2;
  for i=2:numel(Tr_)

    % get centroid of the current cluster
    x = X1_median(idxMaxIn(i),1);
    z = X1_median(idxMaxIn(i),3);

    % compute distances to all accepted clusters
    d_xz = sqrt( (XZ(:,1)-x).^2 + (XZ(:,2)-z).^2 );

    if sum( d_xz<threshNonmax ) > 0, continue; end

    % add current cluster to output struct
    Tr{iTr} = Tr_{idxMaxIn(i)};
    XZ(iTr,:) = [ X1_median(idxMaxIn(i),1), X1_median(idxMaxIn(i),3) ];
    uv1_max{iTr} = uv1_cluster{idxMaxIn(i)};
    uv2_max{iTr} = uv2_cluster{idxMaxIn(i)};

    objects{iTr}.Tr = Tr_{idxMaxIn(i)};
    objects{iTr}.uv = uv1_cluster{idxMaxIn(i)};
    objects{iTr}.X1 = X1_cluster{idxMaxIn(i)};
    objects{iTr}.X2 = X2_cluster{idxMaxIn(i)};
    objects{iTr}.X_median = X1_median(idxMaxIn(i),:);

    iTr = iTr+1;
  end
   
  time_spent = time_spent + toc(t2);
  %fprintf('Object proposals computed in\t%6.2f s\n',time_spent);
    
end