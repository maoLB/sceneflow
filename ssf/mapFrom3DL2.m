function [map,objectLabels,superpixelObjects,tSpent] = ...
  mapFrom3DL2(S,L,objects,Pd,thresh,I,adjMat,Tr,roadMask)
  % Initialization of ks.

  % Sets all superpixels inside the 'roadMask' to background.
  % Computes 3D superpixel centers from current disp planes (not returned).

  % Returns:
  % map                 A vector of initial ks for the superpixels
  % objectLabels        (h, w, nO) matrix containing object masks
  % superpixelObjects   (nS, nO) matrix with binary indicators of possible
  %                     associations between each superpixel and each object

  tStart = tic;

  nO = numel(objects);
  nS = numel(S);

  % initialization of containers
  map = ones(nS,1);

  [h,w] = size(L);
  objectLabels = zeros(h,w,nO);
  superpixelObjects = false(nS,nO+1);
  superpixelObjects(:,1) = 1;
  alphaMap = zeros(h,w);

  % get 3D superpixel centers from current disp planes
  C = superpixelCenters(S,Pd);

  % get road indicator
  ind  = sub2ind(size(L),round([S.cv]'),round([S.cu]'));
  road = roadMask(ind)>0.5;
  C(road,:) = 1000;

  for iO=1:nO

    % get median object coordinates
    Xmedian = median(objects{iO}.X1);

    % compute 3D distance to superpixel centers
    dX   = bsxfun(@minus,C,Xmedian);
    dist = sqrt(dX(:,1).^2+dX(:,2).^2+dX(:,3).^2);

    % get indices of superpixels closer than thresh
    idc  = find(dist<thresh);

    o = zeros(size(L));

    % (first) add all adjacent superpixels
    for i=1:numel(idc)   
      adj = [find(adjMat(idc(i),:))';find(adjMat(:,idc(i)))]; % adjMat is upper tri
      for j=1:numel(adj)
        o(L==adj(j)) = iO;
      end
    end

    % (second) mark the original superpixels
    for i=1:numel(idc)   
      o(L==idc(i)) = iO;
    end

    % apply closing
    se = strel('square',100);
    m  = imclose(o,se) > 0.5; % boolean mask

    o(m) = iO;
    alphaMap(m) = alphaMap(m)+1;

    objectLabels(:,:,iO)=o;
  end

  % set all road superpixels to background
  for iS = 1:numel(S)
    if road(iS)
      for iO = 1:nO
        ind3d = sub2ind(size(objectLabels),S(iS).v,S(iS).u,iO*ones(S(iS).nPix,1));
        objectLabels(ind3d) = 0;
      end
    end
  end

  idx = sub2ind(size(L),round([S.cv]),round([S.cu]));
  for iO = 1:nO
    oL = objectLabels(:,:,iO);
    map(oL(idx)>0) = iO+1;
    superpixelObjects(oL(idx)>0,iO+1) = 1;
  end

  tSpent = toc(tStart);

end
