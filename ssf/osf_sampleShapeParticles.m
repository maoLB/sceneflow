function [shapeParticles,objMap] = ...
  osf_sampleShapeParticles(S,nO,shape,stdShape,shapeParticles,probDist,objLabels,superpixelObjects)
  %SAMPLESHAPEPARTICLES generates particles starting from current shape estimate

  nS  = numel(S);                     % number of superpixels
  [~,~,nP] = size(shapeParticles);    % number of particles
  nPh = round(nP/2);                  % half the number of particles

  [h,w,~]  = size(objLabels);

  objMap = zeros(h,w,nO);             % for visualization

  for iS = 1:nS % all superpixels

    % first particle is current MAP solution
    shapeParticles(iS,:,1) = shape(iS,:,1);

    % Gaußian sampling around current map parameters 
    % for half the number of particles
    for d = 1:3 % alpha, beta, gamma
      shapeParticles(iS,d,2:nPh) = normrnd( shape(iS,d,1), stdShape(d), nPh-1, 1);
    end

    % Propose parameters of nearby superpixels
    idx  = discretesample(probDist(iS,:),nP-nPh);

    for iP = 1:nP-nPh
      abc   = abcFromAlphaBetaGamma( shape( idx(iP) ,1:3 ,1), [S(idx(iP)).cu,S(idx(iP)).cv] );
      alpha = alphaBetaGammaFromAbc( abc, [S(iS).cu,S(iS).cv] );
      shapeParticles(iS,1:3,nPh+iP) = alpha;
    end

    % Sampling of k (object proposal)
    candidates = find(superpixelObjects(iS,:));
    idxO = randi(numel(candidates),[1,nP-1]);

    shapeParticles(iS,4,2:end) = candidates(idxO); 
  end   
    
end
