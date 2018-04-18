function R = computeRoadPlaneTransformation(D,intrinsics)

% get image size
width  = size(D,2);
height = size(D,1);

% road plane segmentation
ransac_roi = [200 width-200 round(height/2) height];
plane_dsi  = computeRoadPlaneEstimateMex(D',int32(ransac_roi),int32(300),single(1));

if plane_dsi(1) < -0.02 || plane_dsi(1) > 0.02 || ...
   plane_dsi(2) < 0.3 || plane_dsi(2) > 0.34 || ...
   plane_dsi(3) < -70 || plane_dsi(3) > -50
    R = false(size(D));
    return;
end

R = segmentRoadPlaneDsi(D,plane_dsi,intrinsics,4);
