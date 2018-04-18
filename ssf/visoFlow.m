function [p_matched,time_spent] = visoFlow( I )
%VISOFLOW Sparse optical flow

  tStart = tic;

  % matching parameters
  param.nms_n                  = 2;   % non-max-suppression: min. distance between maxima (in pixels)
  param.nms_tau                = 50;  % non-max-suppression: interest point peakiness threshold
  param.match_binsize          = 50;  % matching bin width/height (affects efficiency only)
  param.match_radius           = 200; % matching radius (du/dv in pixels)
  param.outlier_flow_tolerance = 5;   % outlier removal: flow tolerance (in pixels)
  param.multi_stage            = 1;   % 0=disabled,1=multistage matching (denser and faster)
  param.half_resolution        = 0;   % 0=disabled,1=match at half resolution, refine at full resolution
  param.refinement             = 0;   % refinement (0=none,1=pixel,2=subpixel)

  param.max_features = 8;
  param.ransac_iters = 40000;

  % init matcher
  matcherMex('init',param);

  % push back images
  matcherMex('push',I{1});
  matcherMex('push',I{2});

  % match features
  matcherMex('match',0);
  p_matched = matcherMex('get_matches',0);

  % close matcher
  matcherMex('close');

  time_spent = toc(tStart);
  % fprintf('Sparse flow computed in\t\t%6.2f s\n',time_spent);

end

