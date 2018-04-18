function [D_gt_occ, D_gt_noc, F_gt_occ, F_gt_noc] = loadOSFGroundTruth(pathDataset,matlabShot)
% Reads ground truth disparity and flow maps

D_gt_occ = cell(1,2);
D_gt_noc = cell(1,2);

for i=0:1

    fn          = sprintf('%06d_10',matlabShot);

    fn_dgt        = fullfile(pathDataset,sprintf('disp_occ_%d/%s.png',i,fn));
    D_gt_occ{i+1} = disp_read(fn_dgt);
    
    fn_dgt        = fullfile(pathDataset,sprintf('disp_noc_%d/%s.png',i,fn));
    D_gt_noc{i+1} = disp_read(fn_dgt);

end

fn_fgt    = fullfile(pathDataset,sprintf('flow_occ/%06d_10.png', matlabShot));
F_gt_occ  = flow_read(fn_fgt);

fn_fgt    = fullfile(pathDataset,sprintf('flow_noc/%06d_10.png', matlabShot));
F_gt_noc  = flow_read(fn_fgt);