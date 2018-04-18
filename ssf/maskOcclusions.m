function [mask] = maskOcclusions(D_in,Pd,Tr)
% project disparity map to object space (using SGM disp), 
% transform 3D points with Tr from ego-motion,
% mask points leaving the image

  % cast disp map for interpolation
  D = single(D_in);

  % get image size
  [h,w] = size(D);

  % interpolate disparity map so all pixels can be projected
  backgroundInterpolationMex(D,D>0);

  % image coordinates in reference frame
  [u,v] = meshgrid(1:w,1:h);

  % project disp map to 3D
  X1 = project([u(:)-1,v(:)-1,D(:)],inv(Pd));

  % transform 3D points according to ego-motion
  X2 = project(X1,Tr);

  % re-project transformed points to left image
  uv2_l = project(X2,Pd);

  % find points leaving the image
  out = uv2_l(:,1)<0 | uv2_l(:,1)>w | uv2_l(:,2)<0 | uv2_l(:,2)>h;

  % mark pixels leaving the image
  iO        = sub2ind(size(D),v(out),u(out));
  mask      = false(size(D));
  mask(iO)  = true;

end