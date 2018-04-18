clear all; close all; dbstop error;
disp('===========================');

% variables: vector specifying the number of states per variable (2=binary)
% factors: cellarray specifying unary, pairwise and high-order potentials
% - v: variables involved in this factor (1-based index)
% - e: energy value for each state combination of the involved variables.
%      For pairwise potentials an energy matrix E can be converted to e
%      via e = E(:)', where the rows of E correspond to the first variable
%      and the columns of E to the second variable. In general, the first
%      entries of e correspond to incrementing the state of the first
%      variable by 1. Here is an example for two variables with 3 states:
%      x1:  1 2 3 1 2 3 1 2 3
%      x2:  1 1 1 2 2 2 3 3 3
%      idx: 1 2 3 4 5 6 7 8 9

% number of states per variable
variables = [3 3 3];

% unary potential
factors{1}.v = 1;
factors{1}.e = [1 -0.5 0];

% pairwise potentials
E = [0 1 1; 1 0 1; 1 1 0];
factors{2}.v = [1 2];
factors{2}.e = E(:)';
factors{3}.v = [2 3];
factors{3}.e = E(:)';

% options
[labels, energy, best_lb] = trwsMex(variables,factors);
labels
energy