%% Enterence of main program 

clear all; close all; clc;  

tether_cfg.problem = 'deployment_problem';
tether_cfg.model_params = {'orientation=vert maxtau=[3,0]'};
tether_cfg.gamma = .98;      % discount factor

if exist('deployment_fuzzyqidemo.mat', 'file'), delete('deployment_fuzzyqidemo.mat'); end;

% Configure the algorithm
cfg = tether_cfg;                      % start from the problem config
cfg.run = 1;                        % set the "run" flag
cfg.datadir = pwd;                  % save the data in the current directory
cfg.datafile = 'deployment_fuzzyqi';  % in this file
cfg.xgrids = symloggrid([5 9 5 7], [0.98 1.2 1.2 1.2]);
cfg.ugrids = symequidgrid([5 1], [3 0]);
cfg.maxiter = 6000;                  % run at most this number of iterations
cfg.eps = 0.01;                      % stop when the difference between consecutive parameter             

fuzzyqi(cfg);

% deployment_plot('fuzzyh datafile=deployment_fuzzyqi fzstyle=qi slice=[NaN,0,NaN,0] figsize=[500,700]');

fuzzyqi('replay datafile=deployment_fuzzyqi x0=[-0.98,0.69,0.116,0] tend=20');



