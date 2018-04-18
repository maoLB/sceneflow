function R = segmentRoadPlaneDsi(D,plane_dsi,intrinsics,threshold)

% get image dimensions
width  = size(D,2);
height = size(D,1);

% horizon
v_min_1 = (-plane_dsi(1)*1-plane_dsi(3))/plane_dsi(2);
v_min_2 = (-plane_dsi(1)*width-plane_dsi(3))/plane_dsi(2);
v_min   = ceil(max(v_min_1,v_min_2))+5;

[U,V] = meshgrid(1:width,1:height);
R = U.*plane_dsi(1)+V.*plane_dsi(2)+plane_dsi(3);
R = single(abs(D-R)<threshold);
R(D<=0) = -1;
R(1:v_min-1,:) = 0;

if ~all(size(R)==size(D))
    R = false(size(D));
    return;
end
