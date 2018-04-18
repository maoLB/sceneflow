% This script demonstrates the usage of osf.m to estimate 3D scene flow as
% described in 'Object Scene Flow for Autonomous Vehicles'

clear all; close all; clc

%% set paths to sourcecode directories
setup_paths

%% set paths to data
pathDataset  = './data';
pathResults  = './results/training';

fnParameters = './ssf/parameters/sf_fast_parameters.mat';

mode = 'training';

[F,D]=osf(pathDataset,pathResults,fnParameters,mode,123);


