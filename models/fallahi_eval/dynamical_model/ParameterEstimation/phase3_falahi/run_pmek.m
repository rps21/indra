% this script will start-up a parallel tempering job for BACCAM3
%
%   to execute from command line on remote server:
%   > nohup matlab -r run_pt > run_pt.out 2> run_pt.err < /dev/null &
  
% path to parallel tempering scripts
% addpath('/home/bobby/Dropbox/MATLAB/Tcell/Parallel_Tempering');
% 
% % path to supplementary distributions
% addpath('/home/bobby/Dropbox/MATLAB/Tcell/Parallel_Tempering/distr/');
% 
% % path to analysis tools (optional)
% addpath('/home/bobby/Dropbox/MATLAB/Tcell/Parallel_Tempering/analysis/');
% 
% addpath('/home/bobby/Dropbox/MATLAB/Tcell/Parallel_Tempering/vis/');
% 
% % path to model-specific files (if not the current directory)
% addpath('/home/bobby/Dropbox/MATLAB/Tcell/Parallel_Tempering/Tcell_model/');
% 
% % path to data (if not the current directory)
% addpath('/home/bobby/Dropbox/MATLAB/Tcell/Parallel_Tempering/Tcell_model/data/');

%load configuration file
cfg = config_pmek();

% start parallel tempering
parallel_tempering(cfg);

%quit;




