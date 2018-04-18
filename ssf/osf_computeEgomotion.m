function [Tr,tSpent] = osf_computeEgomotion( Il, Ir, Pd)
% Computes egomotion from subsequent (grayscale) stereo pairs
% Requires libviso2 to be found in the matlab path

  tStart = tic;
  
  %% matching parameters
  param.nms_n                  = 2;   % non-max-suppression: min. distance between maxima (in pixels)
  param.nms_tau                = 50;  % non-max-suppression: interest point peakiness threshold
  param.match_binsize          = 50;  % matching bin width/height (affects efficiency only)
  param.match_radius           = 200; % matching radius (du/dv in pixels)
  param.match_disp_tolerance   = 1;   % du tolerance for stereo matches (in pixels)
  param.outlier_disp_tolerance = 5;   % outlier removal: disparity tolerance (in pixels)
  param.outlier_flow_tolerance = 5;   % outlier removal: flow tolerance (in pixels)
  param.multi_stage            = 1;   % 0=disabled,1=multistage matching (denser and faster)
  param.half_resolution        = 0;   % 0=disabled,1=match at half resolution, refine at full resolution
  param.refinement             = 1;   % refinement (0=none,1=pixel,2=subpixel)

  %% odometry parameters
  param.max_features = 4;
  param.bucket_width = 50;
  param.bucket_height = 50;
  param.f     = Pd(1,1);
  param.cu    = Pd(1,3);
  param.cv    = Pd(2,3);
  param.base  = Pd(3,4)/Pd(1,1);
  param.ransac_iters = 40000;
  param.inlier_threshold = 3.0;
  param.reweighting = 1;
  
  %% estimate motion
  visualOdometryStereoMex('init',param);
  [~]= visualOdometryStereoMex('process',Il{1},Ir{1});
  Tr = visualOdometryStereoMex('process',Il{2},Ir{2});
  visualOdometryStereoMex('close');
 
  tSpent = toc(tStart);

end
  