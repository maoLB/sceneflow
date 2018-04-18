function [L, D_slic, D_sgm_l, D_sgm_r, B, adjMat, sup_boundaries, Il, Ir] ...
     = loadOSFInput(pathDataset,pathResults,matlabShot)

% Returns:
Il = cell(1,2);       % 2 subsequent frames from left camera
Ir = cell(1,2);       % 2 subsequent frames from right camera
D_slic  = cell(1,2);  % 2 subsequent disparity maps from stereoSLIC (left cam)
D_sgm_l = cell(1,2);  % 2 subsequent disparity maps from SGM (left cam)
D_sgm_r = cell(1,2);  % 2 subsequent disparity maps from SGM (right cam)
L  = cell(1,2);       % 2 subsequent label matrices indexing superpixels
B  = cell(1,2);       % superpixel boundaries from stereoSLIC

for i=0:1

    fn          = sprintf('%06d_1%d',matlabShot,i);

    fn_imgL     = fullfile(pathDataset,sprintf('image_2/%s.png', fn));
    Il{i+1}     = imread(fn_imgL);
    
    fn_imgR     = fullfile(pathDataset,sprintf('image_3/%s.png', fn));
    Ir{i+1}     = imread(fn_imgR);
    
    if 0
      % read precomputed results of SGM and stereoSLIC
      fn_disp     = fullfile(pathResults,sprintf('sgm/%s_left_disparity.png', fn));
      D_sgm_l{i+1}= disp_read(fn_disp);

      fn_disp     = fullfile(pathResults,sprintf('sgm/%s_right_disparity.png', fn));
      D_sgm_r{i+1}= disp_read(fn_disp);
      

      fn_disp     = fullfile(pathResults,sprintf('stereoSLIC/%s_disparity.png', fn));
      D_slic{i+1} = disp_read(fn_disp);
      
      fn_L        = fullfile(pathResults,sprintf('stereoSLIC/%s.png', fn));
      L{i+1}      = imread(fn_L)+1; % correct labels to 1-based indexing

      fn_bnd      = fullfile(pathResults,sprintf('stereoSLIC/%s_boundary.png', fn));
      B{i+1}      = imread(fn_bnd);
    
    end
    
end

% online computation
[D_slic{1}, L{1}] = spsStereoMex(Il{1},Ir{1});
[D_slic{2}, L{2}] = spsStereoMex(Il{2},Ir{2});

D_sgm_l=D_slic;

L{1}=L{1}+1; B=0;
L = uniqueLabels(L);

% compute adjacency matrix
adjMat = superpixel_adjacencies(L{1});

% find boundary pixels
[M,~] = superpixel_masks(L{1});
sup_boundaries = superpixel_boundaries(M, L{1});

end
