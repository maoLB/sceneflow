% Main script calling object scene flow code on the complete 
% test data set. Please download the data from:
% http://www.cvlibs.net/download.php?file=data_scene_flow.zip

clear all; close all; clc;

%% set paths to sourcecode directories
setup_paths

%% set paths to data
pathDataset  = './data';
pathResults  = './results/testing';

fnParameters = './ssf/parameters/sf_fast_parameters.mat'; % fast variant
%fnParameters = './ssf/parameters/sf_parameters.mat'; 		% full method

mode = 'testing';

for i = 0:199
  
  [F,D]=osf(pathDataset,pathResults,fnParameters,mode,i);
  
end
