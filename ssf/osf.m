function [F,D_est]=osf(pathDataset,pathResults,fnParameters,mode,i_shot)

  %% convert shot index to double if necessary
  if ischar(i_shot)
    i_shot = str2double(i_shot);
  end

  %% load parameters
  load(fnParameters);
  p.mode = mode;

  % adjust paths wrt. processing mode
  pathDataset = fullfile(pathDataset,p.mode);

  %% redirect output
  fLog = 1;     % output to cmd window (1)
    
  %% read calibration
  filename = fullfile(pathDataset,'calib',sprintf('%06d.txt',i_shot));
  [ P0,P1,K,Pd,m,b_m ] = loadStereoCalib(filename);

  %% read input data
  [L,D_slic,D_sgm,D_sgm_r,~,adjMat,sup_boundaries,Il,Ir] = loadOSFInput(pathDataset,pathDataset,i_shot);
  [h,w,~] = size(Il{1});

  %% in training mode read ground truth
%   if strcmp(p.mode,'training')
%     [D_gt_occ, D_gt_noc, F_gt_occ, F_gt_noc] = loadOSFGroundTruth(pathDataset,i_shot);
%   end
  
  %% init random number generator for deterministic results
  s = RandStream('mt19937ar','Seed',0);
  RandStream.setGlobalStream(s);
  
  %% convert label matrix to segmentation struct
  [S,Dp,locDist] = labellingToStruct(L{1},D_slic{1});
  probDist       = locDistToProb(locDist,100);
  nS = numel(S);
  
  %% compute grayscale images
  [Il_gray,Ir_gray,tGray] = grayImages(Il,Ir);
  
  %% compute census transformed images
  [Cl,Cr,tCens] = censusTransformAll(Il_gray,Ir_gray);
  
  %% estimate roadplane from disparity map
  roadMask  = computeRoadPlaneTransformation(single(D_slic{1}),[P0(1,1),P0(2,2),P0(1,3),P0(2,3),b_m]);

  % erode road mask
  se = strel('disk',10,0);
  roadMask = imerode(roadMask,se);
  
  %% compute sparse optical flow
  [sparseFlowLeft,tFlow]   = visoFlow(Il_gray);
  sparseFlowLeft_all       = sparseFlowLeft;
  F_sparse{1} = matchesToFlowMap(sparseFlowLeft,w,h);
  
  %% estimate ego-motion
  [Tr_fldp,tMotion] = osf_computeEgomotion(Il_gray, Ir_gray, Pd);
  
  %% mask points leaving the image
  mask_l = maskOcclusions(D_sgm{1},Pd,Tr_fldp);
  mask3D_l = cat(3,mask_l,mask_l,mask_l);
  F_sparse{1}(mask3D_l) = 0;
  
  %% cross-term observations (connecting Il{1} and Ir{2})
  % from combination of flow & disp (in that order)
  sparseCrossTerm2   = crossTerm_flowDisp(S,D_sgm{2},F_sparse{1});
  F_C1_flowDisp{1}   = matchesToFlowMap(sparseCrossTerm2,w,h);      
  
  %% get object hypotheses
  shape(:,:,1)    = [S.alpha; S.beta; S.gamma; ones(1,nS)]';

  euler = EulerFromRot(Tr_fldp(1:3,1:3));
  object(1,:,1) = [euler, Tr_fldp(1:3,4)'];

  [F,motionMap] = flowMapFromShapeAndMotion(S,shape(:,:,1),object(:,:,1),h,w,K,b_m);
  
  [Tr,objects,tOP] = objectProposalsFlowDisp(sparseFlowLeft_all,D_sgm,Il,F,Pd,s);
  
  %% ==================== Inference ========================
  
  %% set up the model
  % variables: vector specifying the number of states per variable (2=binary)
  % variables = [ superpixels, objects ]
  variables = [ ones(1, nS) * p.nShapeParticles, ones(1, p.nO) * p.nObjectParticles ];
  nV = length(variables); % number of variables
  
  % shape: contains final shape parameters for each superpixel
  % first index refers to the superpixel
  % second index references alpha, beta, gamma (in this order) and the object indicator k
  % third index corresponds to the number of iterations starting with
  % the initialization
  shape           = zeros(nS, 4, p.nIters+1);
  shape(:,:,1)    = [S.alpha; S.beta; S.gamma; ones(1,nS)]';

  % shapeParticles contains the current proposal planes
  % starting with the previous MAP solution
  shapeParticles  = zeros(nS, 4, p.nShapeParticles);
  
  % object: contains final motion parameters for each object
  % first index refers to the object
  % second index references rotations and translations in x,y,z (in this order)
  % third index corresponds to the number of iterations starting with
  % the initialization
  object          = zeros(p.nO, 6, p.nIters+1);
  
  % background motion
  euler = EulerFromRot(Tr_fldp(1:3,1:3));
  object(1,:,1) = [euler, Tr_fldp(1:3,4)'];

  % objectParticles contains the current motion proposals
  % starting with the previous MAP solution
  objectParticles = zeros(p.nO, 6, p.nObjectParticles);

  % containers for optimization results
  energy  = zeros(p.nIters+1,1);
  tIter   = zeros(p.nIters+1,1);
  map     = zeros(p.nIters+1,nV);
  
  % init objects if there are hypotheses
  if numel(Tr) > 0
    [mapOfObj,objLabels,superpixelObjects,tInit] = ...
      mapFrom3DL2(S,L{1},objects(1:min(p.nO-1,numel(Tr))),Pd,2.5,Il{1},adjMat,Tr,roadMask);
    
    % init shape particles with a plausible k
    shape(:,4,1) = mapOfObj;

    % init the chosen number of objects
    for iOP = 1:min(p.nO-1,numel(Tr))
      euler = EulerFromRot(Tr{iOP}(1:3,1:3));
      object(iOP+1,:,1) = [euler, Tr{iOP}(1:3,4)'];
    end
    
  else % there are no hypotheses
    mapOfObj  = ones(nS,1);
    objLabels = zeros(h,w);
    superpixelObjects = ones(nS,1);
  end

  %% optimization loop
  for i = 2:p.nIters+1

    % start timer
    tI = tic;

    % reset particles
    shapeParticles  = zeros(size(shapeParticles));
    objectParticles = zeros(size(objectParticles));

    % sample particles
    tParticles = tic;
    
    [shapeParticles,objMap] = osf_sampleShapeParticles(S,p.nO,shape(:,:,i-1),p.stdShape(i-1,:),shapeParticles,probDist,objLabels,superpixelObjects);
    [objectParticles,~]     = osf_sampleObjectParticles(object(:,:,i-1),p.stdObject(i-1,:),objectParticles,p.stdBackground);
    
    %% pairwise data term
    
    % init pairwise factors
    factors_data = cell(1,nS*p.nO);
    for iS = 1:nS % all superpixels
      for iO = 1:p.nO % all objects
        factors_data{(iS-1)*p.nO+iO}.v = [iS,nS+iO];
        factors_data{(iS-1)*p.nO+iO}.e = zeros(p.nShapeParticles,p.nObjectParticles);
      end
    end
    
    % compute relevant data terms
    tData = tic;
    E_data = computeDataTermMex(S,K,m,shapeParticles,objectParticles,p,D_sgm{1},F_sparse{1},F_C1_flowDisp{1},Cl,Cr,0);

    for iS = 1:nS % all superpixels

      for iSP = 1:p.nShapeParticles

        k  = shapeParticles(iS,4,iSP);

        if k==1
          factors_data{(iS-1)*p.nO+k}.e(iSP,:) = E_data( iS, (iSP-1)*p.nObjectParticles+1:(iSP-1)*p.nObjectParticles+p.nObjectParticles )-p.bgPrior;
        else
          factors_data{(iS-1)*p.nO+k}.e(iSP,:) = E_data( iS, (iSP-1)*p.nObjectParticles+1:(iSP-1)*p.nObjectParticles+p.nObjectParticles );
        end     
      end
    end
    
    %% reshape factors to row vectors
    for iS = 1:nS % all superpixels
      for iO = 1:p.nO % all objects
        factors_data{(iS-1)*p.nO+iO}.e = factors_data{(iS-1)*p.nO+iO}.e(:)';
      end
    end
    
    %% smoothness term
    [E_pw_disp,E_pw_norm,E_pw_pott,pairs] = computeSmoothnessTermMex(S,adjMat,sup_boundaries,shapeParticles,p.nShapeParticles,h,K,m,p,0);
    factors_pw = pairwiseFactors(pairs,E_pw_disp,p.wDispSmooth,E_pw_norm,p.wNormalSmooth,E_pw_pott,p.wPottsSmooth);
    
    %% optimization and update
    factors = [factors_data,factors_pw];
    
    options.maxIter = 200;
    [map(i,:), energy(i)] = trwsMex(variables, factors, options);

    % save current MAP shape
    for iS = 1:nS
      shape(iS,:,i)  = shapeParticles(iS,:,map(i,iS));
    end

    % save current MAP objects
    for iV = nS+1:nS+p.nO
      iO = iV-nS; % 1-based object index
      object(iO,:,i) = objectParticles(iO,:,map(i,iV));
    end
    
    % update std deviation of particle distribution
    p.stdShape(i,:)   = [ p.stdShape(1,1)  * exp(-(i-1)/p.decayShape), p.stdShape(1,2)   * exp(-(i-1)/p.decayShape), p.stdShape(1,3)   * exp(-(i-1)/p.decayShape) ];
    p.stdObject(i,:)  = [ p.stdObject(1,1) * exp(-(i-1)/p.decayMotion), p.stdObject(1,2) * exp(-(i-1)/p.decayMotion), p.stdObject(1,3) * exp(-(i-1)/p.decayMotion),...
                          p.stdObject(1,4) * exp(-(i-1)/p.decayMotion), p.stdObject(1,5) * exp(-(i-1)/p.decayMotion), p.stdObject(1,6) * exp(-(i-1)/p.decayMotion) ];
    p.stdBackground(i,:) = p.stdBackground(1,:).* exp(-(i-1)/p.decayMotion);
    
    tIter(i) = toc(tI);
    fprintf(fLog,'Iteration %2d finished in % 6.2f s',i-1,tIter(i));
    fprintf(fLog,', final energy: % 10.2f\n',energy(i));
    
    % disparity map from current shape parameters
    D1 = dispMapFromPlanes(S,shape(:,1:3,i),h,w);
    [dD,~] = computeDisparityFlow(S,shape(:,:,i),object(:,:,i),h,w,K,b_m );

    % flow map from current shape and motion parameters
    [F,motionMap] = flowMapFromShapeAndMotion(S,shape(:,:,i),object(:,:,i),h,w,K,b_m);

    D_est{1} = D1;
    D_est{2} = D1+dD;
    
    D_est{1}(D_est{1}<0)=0;
    D_est{2}(D_est{2}<0)=0;
    
  %[nBad_occ,nValid_occ,pct_occ_iter(i-1),idx_occ] = getBadPixels(D_est,F,D_gt_occ,F_gt_occ,p.tau_px,p.tau_pct);
    
  end % optimization loop
  
  %% save results
  if ~exist( fullfile(pathResults,'disp_0'), 'dir'), mkdir(pathResults,'disp_0'); end;
  if ~exist( fullfile(pathResults,'disp_1'), 'dir'), mkdir(pathResults,'disp_1'); end;
  if ~exist( fullfile(pathResults,'flow'),   'dir'), mkdir(pathResults,'flow');   end;
  if ~exist( fullfile(pathResults,'label'),  'dir'), mkdir(pathResults,'label');  end;

  filenameDisp = fullfile(pathResults,sprintf('disp_0/%06d_10.png',i_shot));
  disp_write(D_est{1},filenameDisp);

  filenameDisp = fullfile(pathResults,sprintf('disp_1/%06d_10.png',i_shot));
  disp_write(D_est{2},filenameDisp);

  filenameFlow = fullfile(pathResults,sprintf('flow/%06d_10.png',i_shot));
  flow_write(F,filenameFlow);

  filenameMap  = fullfile(pathResults,sprintf('label/%06d_10.png',i_shot));
  imwrite(uint8(motionMap),filenameMap);
  
end
  
  