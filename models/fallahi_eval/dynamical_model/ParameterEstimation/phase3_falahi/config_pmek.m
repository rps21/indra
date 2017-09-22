function [cfg] = config_bigmech_test()
% CONFIG_BACCAM3 configure BACCAM3 viral model for parallel tempering
%

% parameter definitions
%   Parameter Priors
%   --------------------------------------------------------------
%   'point'       : argument is 'value'
%   'uniform'     : arguments are 'min' and 'max' (must be finite)
%   'laplace'     : arguments are 'mu' (mean) and 'b' (width)
%   'normal'      : arguments are 'mu' and 'sigma'
%   'lognormal'   : arguments are 'mu' and 'sigma'
%   'beta'        : arguments are 'a' and 'b'
%   'exponential' : argument is 'beta' (mean)
% 


cfg.param_defs = { ...

  struct('name','grb2_lox','prior','point','value',log10(10000*8.7223041193),'units','?'), ... 
  struct('name','mapk1_lox','prior','point','value',log10(10000*5.8336394997),'units','?'), ... 
  struct('name','kras_lox','prior','point','value',log10(10000*6.0946099346),'units','?'), ... 
  struct('name','sos1_lox','prior','point','value',log10(10000*8.1375909542),'units','?'), ... 
  struct('name','map2k1_lox','prior','point','value',log10(10000*8.1375909542),'units','?'), ... 
  struct('name','egfr_lox','prior','point','value',log10(10000*4.7220864391),'units','?'), ...
  struct('name','mapk3_lox','prior','point','value',log10(10000*6.3067137199),'units','?'), ... 
  struct('name','braf_lox','prior','point','value',log10(10000*4.3408439522),'units','?'), ... 
  struct('name','egf_lox','prior','point','value',log10(1000*2.7983094633),'units','?'), ... 
  struct('name','phos_lox','prior','uniform','min',log10(1e2),'max',log10(1e6),'units','?'), ... 
  struct('name','kf_kb_act_1','prior','uniform','min',log10(1e-7),'max',log10(1e-5),'units','?'), ... 
  struct('name','kf_e_autophos_1','prior','uniform','min',log10(1e-1),'max',log10(11e1),'units','?'), ... 
  struct('name','kf_eg_bind_1','prior','uniform','min',log10(1e-7),'max',log10(5e-5),'units','?'), ... 
  struct('name','kr_eg_bind_1','prior','uniform','min',log10(1e-3),'max',log10(1e-1),'units','?'), ... 
  struct('name','kf_ee_bind_1','prior','uniform','min',log10(1e-7),'max',log10(1e-5),'units','?'), ... 
  struct('name','kr_ee_bind_1','prior','uniform','min',log10(1e-3),'max',log10(1e-1),'units','?'), ... 
  struct('name','kf_gs_bind_1','prior','uniform','min',log10(1e-7),'max',log10(1e-5),'units','?'), ... 
  struct('name','kr_gs_bind_1','prior','uniform','min',log10(1e-3),'max',log10(1e-1),'units','?'), ... 
  struct('name','kf_sk_gef_1','prior','uniform','min',log10(1e-7),'max',log10(1e-5),'units','?'), ... 
  struct('name','kf_bm_bind_1','prior','uniform','min',log10(1e-7),'max',log10(1e-5),'units','?'), ... 
  struct('name','kr_bm_bind_1','prior','uniform','min',log10(1e-3),'max',log10(1e-1),'units','?'), ... 
  struct('name','kc_bm_phosphorylation_1','prior','uniform','min',log10(1e-1),'max',log10(1e1),'units','?'), ... 
  struct('name','kf_bm_bind_2','prior','uniform','min',log10(1e-7),'max',log10(1e-5),'units','?'), ... 
  struct('name','kr_bm_bind_2','prior','uniform','min',log10(1e-3),'max',log10(1e-1),'units','?'), ... 
  struct('name','kc_bm_phosphorylation_2','prior','uniform','min',log10(1e-1),'max',log10(1e1),'units','?'), ... 
  struct('name','kf_bm_bind_3','prior','uniform','min',log10(1e-7),'max',log10(1e-5),'units','?'), ... 
  struct('name','kr_bm_bind_3','prior','uniform','min',log10(1e-3),'max',log10(1e-1),'units','?'), ... 
  struct('name','kc_bm_phosphorylation_3','prior','uniform','min',log10(1e-1),'max',log10(1e1),'units','?'), ... 
  struct('name','kf_mm_bind_1','prior','uniform','min',log10(1e-7),'max',log10(1e-5),'units','?'), ... 
  struct('name','kr_mm_bind_1','prior','uniform','min',log10(5e-3),'max',log10(5e-1),'units','?'), ... 
  struct('name','kc_mm_phosphorylation_1','prior','uniform','min',log10(1e-1),'max',log10(1e1),'units','?'), ... 
  struct('name','kf_mm_bind_2','prior','uniform','min',log10(5e-7),'max',log10(5e-5),'units','?'), ... 
  struct('name','kr_mm_bind_2','prior','uniform','min',log10(1e-3),'max',log10(1e-1),'units','?'), ... 
  struct('name','kc_mm_phosphorylation_2','prior','uniform','min',log10(1e-1),'max',log10(1e1),'units','?'), ... 
  struct('name','kf_mb_bind_1','prior','uniform','min',log10(1e-7),'max',log10(1e-5),'units','?'), ... 
  struct('name','kr_mb_bind_1','prior','uniform','min',log10(1e-3),'max',log10(1e-1),'units','?'), ... 
  struct('name','kc_mb_phosphorylation_1','prior','uniform','min',log10(1e-1),'max',log10(1e1),'units','?'), ... 
  struct('name','kf_mm_bind_3','prior','uniform','min',log10(1e-7),'max',log10(1e-5),'units','?'), ... 
  struct('name','kr_mm_bind_3','prior','uniform','min',log10(1e-7),'max',log10(1e-5),'units','?'), ... 
  struct('name','kc_mm_phosphorylation_3','prior','uniform','min',log10(1e-1),'max',log10(1e1),'units','?'), ... 
  struct('name','kf_mm_bind_4','prior','uniform','min',log10(1e-7),'max',log10(1e-5),'units','?'), ... 
  struct('name','kr_mm_bind_4','prior','uniform','min',log10(1e-3),'max',log10(1e-1),'units','?'), ... 
  struct('name','kc_mm_phosphorylation_4','prior','uniform','min',log10(1e-1),'max',log10(1e1),'units','?'), ... 
  struct('name','kf_mb_bind_2','prior','uniform','min',log10(1e-7),'max',log10(1e-5),'units','?'), ... 
  struct('name','kr_mb_bind_2','prior','uniform','min',log10(1e-3),'max',log10(1e-1),'units','?'), ... 
  struct('name','kc_mb_phosphorylation_2','prior','uniform','min',log10(1e-2),'max',log10(1),'units','?'), ... 
  struct('name','kf_mm_bind_5','prior','uniform','min',log10(1e-7),'max',log10(1e-5),'units','?'), ...
  struct('name','kr_mm_bind_5','prior','uniform','min',log10(1e-3),'max',log10(1e-1),'units','?'), ... 
  struct('name','kc_mm_phosphorylation_5','prior','uniform','min',log10(1e-1),'max',log10(1e1),'units','?'), ... 
  struct('name','kf_mm_bind_6','prior','uniform','min',log10(1e-7),'max',log10(1e-5),'units','?'), ... 
  struct('name','kr_mm_bind_6','prior','uniform','min',log10(1e-3),'max',log10(1e-1),'units','?'), ... 
  struct('name','kc_mm_phosphorylation_6','prior','uniform','min',log10(1e-2),'max',log10(1e1),'units','?'), ... 
  struct('name','kf_drug_braf','prior','uniform','min',log10(1e-8),'max',log10(1e-5),'units','?'), ... 
  struct('name','kr_drug_braf','prior','uniform','min',log10(1e-4),'max',log10(1e-1),'units','?'), ... 
  struct('name','mapk1_p_decay_slow','prior','uniform','min',log10(1e-6),'max',log10(1e-3),'units','?'), ... 
  struct('name','kf_erk_phos','prior','uniform','min',log10(1e-6),'max',log10(1e-3),'units','?'), ... 
  struct('name','kr_erk_phos','prior','uniform','min',log10(1e-3),'max',log10(5e-1),'units','?'), ... 
  struct('name','kc_erk_phos','prior','uniform','min',log10(1e-1),'max',log10(1e1),'units','?'), ... 
  struct('name','kf_drug_phos','prior','uniform','min',log10(1e-3),'max',log10(1e2),'units','?'), ... 
  struct('name','kr_drug_phos','prior','uniform','min',log10(1e-8),'max',log10(1e-3),'units','?'), ... 
   struct('name','map2k1_decay_slow','prior','uniform','min',log10(1e-6),'max',log10(1e-2),'units','?'), ...
  struct('name','grb2_rvh','prior','point','value',log10(10000*8.452894354),'units','?'), ... 
  struct('name','mapk1_rvh','prior','point','value',log10(10000*6.5677369431),'units','?'), ... 
  struct('name','kras_rvh','prior','point','value',log10(10000*5.3982172252),'units','?'), ... 
  struct('name','sos1_rvh','prior','point','value',log10(10000*7.0660732056),'units','?'), ... 
  struct('name','map2k1_rvh','prior','point','value',log10(10000*7.0660732056),'units','?'), ... 
  struct('name','egfr_rvh','prior','point','value',log10(10000*3.5447896993),'units','?'), ...
  struct('name','mapk3_rvh','prior','point','value',log10(10000*5.874877607),'units','?'), ... 
  struct('name','braf_rvh','prior','point','value',log10(10000*3.8434139219),'units','?'), ... 
  struct('name','egf_rvh','prior','point','value',log10(10000*3.5447896993),'units','?'), ... 
  struct('name','phos_rvh','prior','uniform','min',log10(1e-2),'max',log10(1),'units','?'), ... 
    

};

% initialize parameter distributions
cfg = init_parameter_defs( cfg.param_defs, cfg );


% observable definitions
cfg.obsv_defs = { ...
  struct('name','perk', 'units','log10'), ...
  struct('name','active_braf', 'units','log10'), ...
  struct('name','free_phosphotase', 'units','log10'), ...
  struct('name','drug_phosphotase', 'units','log10'), ...
    struct('name','drug_total', 'units','log10'), ...
  struct('name','drug_bound', 'units','log10'), ...

};

% initialize observable structs
cfg = init_observable_defs( cfg.obsv_defs, cfg );


% custom model configuration
cfg.sim_tstart = 0;           % simulation start time, time units
cfg.sim_tstop  = 172800;          % simulation stop time, time units (seconds)
cfg.sim_dt     = 10;         % time step for trajectory, time units (seconds)
cfg.time_units = 'seconds';      % time units
cfg.big_energy = 1e29;        % large energy value (arbitrary units
cfg.tolerance  = 1e-4;        % tolerance factor for comparing things to zero
cfg.maxlogH    = log10(1e6);  % maximum value of 'H' at end of simulation
cfg.timepenalty = 10;         % penalty for long integration times

cfg.plot_xlim = [cfg.sim_tstart cfg.sim_tstop];   % x-axis limits for plotting experiments
cfg.plot_ylims = {[-1 9], [-1 9], [-2 7]};        % y-axis limits for plotting experiments


% load experimental data

load('exp_time_pts.mat')
load('pmek_lox_mean.mat')
load('pmek_lox_std.mat')
load('pmek_rhv_mean.mat')
load('pmek_rhv_std.mat')

% load('testdata_weight.mat')



% pten_l = load('pten_lowdose_normalized');
% pten_h = load('pten_highdose_normalized');

% pten_lowdose = pten_l.pten_lowdose_normalized;
% pten_highdose = pten_h.pten_highdose_normalized;


% experiment name
cfg.data{1}.name = 'low dose';
cfg.data{2}.name = 'mid dose';
cfg.data{3}.name = 'high dose';
cfg.data{4}.name = 'low dose';
cfg.data{5}.name = 'mid dose';
cfg.data{6}.name = 'high dose'

% set up experimental timepoints and map

% timepoints = pten_lowdose(1:9,1);
timepoints = exp_time_pts;
cfg.timepoints = timepoints;
cfg.data{1}.time = [cfg.sim_tstart : cfg.sim_dt : cfg.sim_tstop]'; 
cfg.data{2}.time = [cfg.sim_tstart : cfg.sim_dt : cfg.sim_tstop]'; 
cfg.data{3}.time = [cfg.sim_tstart : cfg.sim_dt : cfg.sim_tstop]'; 



cfg.data{1}.time_map = struct(  ...
     't1',    find(cfg.data{1}.time == timepoints(1)),  ...
     't2',    find(cfg.data{1}.time == timepoints(2)),  ...
     't3',    find(cfg.data{1}.time == timepoints(3)),  ...
     't4',    find(cfg.data{1}.time == timepoints(4)),  ...
     't5',    find(cfg.data{1}.time == timepoints(5)) );
%      't6',    find(cfg.data{1}.time == timepoints(6)),  ...
%      't7',    find(cfg.data{1}.time == timepoints(7)),  ...
%      't8',    find(cfg.data{1}.time == timepoints(8)),  ...
%      't9',    find(cfg.data{1}.time == timepoints(9))  ...
% );

cfg.data{2}.time_map = struct(  ...
     't1',    find(cfg.data{1}.time == timepoints(1)),  ...
     't2',    find(cfg.data{1}.time == timepoints(2)),  ...
     't3',    find(cfg.data{1}.time == timepoints(3)),  ...
     't4',    find(cfg.data{1}.time == timepoints(4)),  ...
     't5',    find(cfg.data{1}.time == timepoints(5)) );
%      't6',    find(cfg.data{1}.time == timepoints(6)),  ...
%      't7',    find(cfg.data{1}.time == timepoints(7)),  ...
%      't8',    find(cfg.data{1}.time == timepoints(8)),  ...
%      't9',    find(cfg.data{1}.time == timepoints(9))  ...
% );

cfg.data{3}.time_map = struct(  ...
     't1',    find(cfg.data{1}.time == timepoints(1)),  ...
     't2',    find(cfg.data{1}.time == timepoints(2)),  ...
     't3',    find(cfg.data{1}.time == timepoints(3)),  ...
     't4',    find(cfg.data{1}.time == timepoints(4)),  ...
     't5',    find(cfg.data{1}.time == timepoints(5)) );
%      't6',    find(cfg.data{1}.time == timepoints(6)),  ...
%      't7',    find(cfg.data{1}.time == timepoints(7)),  ...
%      't8',    find(cfg.data{1}.time == timepoints(8)),  ...
%      't9',    find(cfg.data{1}.time == timepoints(9))  ...
% );


% Set dose for each initial condition
cfg.data{1}.dose = 1; 
cfg.data{2}.dose = 2; 
cfg.data{3}.dose = 3; 
cfg.data{4}.dose = 1; 
cfg.data{5}.dose = 2;
cfg.data{6}.dose = 3; 

%%%%%%%%%% Low Dose %%%%%%%%%%%%%

% cfg.data{1}.mean = nan(length(fieldnames(cfg.data{1}.time_map)), cfg.nobsv);
% cfg.data{1}.mean(:,1) = pten_lowdose(1:9,2); 
% cfg.data{1}.mean(:,2) = nanmean(pten_mrna_lowdose(1:9,2:3)'); 
% cfg.data{1}.mean(:,3) = nanmean(akt_p308_lowdose(1:9,2:3)'); 
% cfg.data{1}.mean(:,4) = akt_p473_lowdose(1:9,2); 
%     cfg.data{1}.mean(8,4) = NaN;    %Only one replicate
% cfg.data{1}.mean(:,5) = foxo_phos_lowdose(1:9,2);

cfg.data{1}.mean = pmek_lox_mean(:,1);
cfg.data{2}.mean = pmek_lox_mean(:,4);
cfg.data{3}.mean = pmek_lox_mean(:,6);
cfg.data{4}.mean = pmek_rhv_mean(:,1);
cfg.data{5}.mean = pmek_rhv_mean(:,4);
cfg.data{6}.mean = pmek_rhv_mean(:,7);

%%%%%%%%%%% High Dose %%%%%%%%%%%
% cfg.data{2}.mean = nan(length(fieldnames(cfg.data{1}.time_map)), cfg.nobsv);
% cfg.data{2}.mean(:,1) = pten_highdose(1:9,2); 
% cfg.data{2}.mean(:,2) = nanmean(pten_mrna_highdose(1:9,2:3)'); 
% cfg.data{2}.mean(:,3) = nanmean(akt_p308_highdose(1:9,2:3)'); 
% cfg.data{2}.mean(:,4) = nanmean(akt_p473_highdose(1:9,2:4)'); 
%     cfg.data{2}.mean(8,4) = NaN;    %Only one replicate
% cfg.data{2}.mean(:,5) = foxo_phos_highdose(1:9,2); 
% 



%%%%%%%%%%% Antigen Removal %%%%%%%%%%%
% cfg.data{3}.mean = nan(length(cfg.data{3}.time), cfg.nobsv);


% % experimental data stdev (ntime x nobsv)
%%%%%%%%%%%% Low Dose %%%%%%%%%%%%
% cfg.data{1}.stdev = nan(length(fieldnames(cfg.data{1}.time_map)), cfg.nobsv);
% cfg.data{1}.stdev(:,1) = pten_lowdose(1:9,3);
% cfg.data{1}.stdev(:,2) = nanstd(pten_mrna_lowdose(1:9,2:3)'); %0.2*cfg.data{1}.mean(:,2);
% cfg.data{1}.stdev(:,3) = nanstd(akt_p308_lowdose(1:9,2:3)'); %0.2*cfg.data{1}.mean(:,2);
% cfg.data{1}.stdev(:,4) = nanstd(akt_p473_highdose(1:9,2:4)'); 
% cfg.data{1}.stdev(:,5) = 0.1; %0.2*cfg.data{1}.mean(:,3);
% 
% ind = find(cfg.data{1}.stdev<=0.025);
% cfg.data{1}.stdev(ind) = .025;
% 
% cfg.data{1}.stdev(:,:) = 0.1;

cfg.data{1}.stdev = pmek_lox_std(:,1);
cfg.data{2}.stdev = pmek_lox_std(:,4);
cfg.data{3}.stdev = pmek_lox_std(:,7);
cfg.data{4}.stdev = pmek_rhv_std(:,1);
cfg.data{5}.stdev = pmek_rhv_std(:,4);
cfg.data{6}.stdev = pmek_rhv_std(:,7);


% %%%%%%%%%%%%% High Dose %%%%%%%%%%%%%%%%
% cfg.data{2}.stdev = nan(length(fieldnames(cfg.data{1}.time_map)), cfg.nobsv);
% cfg.data{2}.stdev(:,1) = pten_highdose(1:9,3);
% cfg.data{2}.stdev(:,2) = nanstd(pten_mrna_highdose(1:9,2:3)'); %0.2*cfg.data{1}.mean(:,2);
% cfg.data{2}.stdev(:,3) = nanstd(akt_p308_highdose(1:9,2:3)'); %0.2*cfg.data{2}.mean(:,2);
% cfg.data{2}.stdev(:,4) = 0.1; %0.2*cfg.data{2}.mean(:,3);
% cfg.data{2}.stdev(:,5) = foxo_phos_highdose(1:9,3); %0.2*cfg.data{2}.mean(:,3);
% 
% ind = find(cfg.data{1}.stdev<=0.025);
% cfg.data{2}.stdev(ind) = .025;
% 
% cfg.data{2}.stdev(:,:) = 0.1;
% %%%%%%%%%%%%% Antigen Removal %%%%%%%%%%%%%%%%
% % cfg.data{3}.stdev = nan(length(cfg.data{3}.time), cfg.nobsv);


% % weights for experimental measurements (ntime x nobsv)
%%%%%%%%%%%%% Low Dose %%%%%%%%%%%%%%%%
% cfg.data{1}.weight = ones( length(fieldnames(cfg.data{1}.time_map)), cfg.nobsv );
% 
% cfg.data{1}.weight(:,1:5) = 0; 
% % cfg.data{1}.weight(:,1:5) = 1; 
% % 
% cfg.data{1}.weight(:,1) = 5;
% % cfg.data{1}.weight(:,1) = 1;    %Low Dose PTEN
% % cfg.data{1}.weight(2,1) = 0;    %Low dose spike in PTEN
% cfg.data{1}.weight(3,1) = 15;   %Low Dose PTEN initial dip
% cfg.data{1}.weight(6,1) = 10;    %Low Dose PTEN rebound
% % % cfg.data{1}.weight(7:9,1) = 7.5;
% % % % % % 
% % cfg.data{1}.weight(:,2) = 2;    %Low Dose mrna
% % % cfg.data{1}.weight(6,2) = 10;    %Low Dose mrna middle
% % % cfg.data{1}.weight(8,2) = 5;    %Low Dose mrna min
% % % % % % 
% % cfg.data{1}.weight(:,3) = 0;    %Low Dose Akt 308
% cfg.data{1}.weight(:,3) = 5;    %Low Dose First min - 0
% % % cfg.data{1}.weight(4,3) = 5;    %Low Dose First min - 0
% % % cfg.data{1}.weight(6,3) = 5;    %Low Dose First min - 0
% % % % cfg.data{1}.weight(7,3) = 5;    %Low Dose First min - 0
% % % 
% % % % % % 
% cfg.data{1}.weight(:,4) = 2;    %Low Dose Akt 473
% % % 
% cfg.data{1}.weight(:,5) = 2;    %Low Dose Foxo phosphoryaltion
% % % cfg.data{1}.weight(1:5,5) = 8;
% % % cfg.data{1}.weight(7:8,5) = 2;

cfg.data{1}.weight = ones(5,1);
cfc.data{1}.weight(1) = 5;
cfc.data{1}.weight(2:end) = 0;
% cfc.data{1}.weight(5) = 10;

cfg.data{2}.weight = ones(5,1);
cfg.data{3}.weight = ones(5,1);
cfg.data{4}.weight = ones(5,1);
cfg.data{5}.weight = ones(5,1);
cfg.data{6}.weight = ones(5,1);

%%%%%%%%%%%%% High Dose %%%%%%%%%%%%%%%%
% % weights for experimental measurements (ntime x nobsv)
% cfg.data{2}.weight = ones( length(fieldnames(cfg.data{1}.time_map)), cfg.nobsv );
% 
% cfg.data{2}.weight(:,1:5) = 0;
% % cfg.data{2}.weight(:,1:5) = 1;
% % % 
% cfg.data{2}.weight(:,1) = 5;    %High dose PTEN
% % % cfg.data{2}.weight(2,1) = 5;    %High dose PTEN initial dip
% % % cfg.data{2}.weight(5,1) = .1;   %Small stdev
% % % cfg.data{2}.weight(9,1) = 10;    %High dose PTEN late min
% % % % % % % % % % 
% % cfg.data{2}.weight(:,2) = 1;    %High dose mrna
% % % cfg.data{2}.weight(6,2) = 10;    %High dose mrna
% % % % cfg.data{2}.weight(8,2) = 5;    %High dose mrna min
% % % % % % % % % 
% % cfg.data{2}.weight(:,3) = 0;    %High Dose Akt 308
% cfg.data{2}.weight(:,3) = 2;    %First min
% % % cfg.data{2}.weight(6,3) = 0.5;    %First min
% % % 
% % % % % % % % 
% cfg.data{2}.weight(:,4) = 2;    %High Dose Akt 473
% % % cfg.data{2}.weight(2,4) = 4;    %High Dose peak
% % % cfg.data{2}.weight(3,4) = 2;    %High Dose first min
% % % cfg.data{2}.weight(4,4) = 7;    %High Dose second peak
% % % cfg.data{2}.weight(6,4) = 5;    %High Dose second min 
% % % cfg.data{2}.weight(7,4) = 2;    %High Dose second min rebound
% % % cfg.data{2}.weight(8,4) = 10;    %High Dose second min rebound
% % % cfg.data{2}.weight(9:end,4) = 0;   
% % % 
% cfg.data{2}.weight(:,5) = 2;    %High dose foxo phosphorylation
% % % cfg.data{2}.weight(2,5) = 5;
% % % cfg.data{2}.weight(3,5) = 5;
% % % cfg.data{2}.weight(4,5) = 7.5;
% % % cfg.data{2}.weight(5,5) = 5;
% % % % cfg.data{2}.weight(6,5) = 2;    %2 hour dip
% % % cfg.data{2}.weight(7,5) = 5;
% % % cfg.data{2}.weight(8:end,5) = 0;

% %%%%%%%%%%%%% High Dose %%%%%%%%%%%%%%%%
% cfg.data{3}.weight = ones( length(cfg.data{3}.time), cfg.nobsv );
% 
% % extract default initial conditions from data file
% %Set in Simulate file now
% 
% get number of experiments
cfg.nexpt = length(cfg.data);


% parallel tempering options
cfg.jobname = 'bigmech_perk';            % job name, for file input/output
cfg.shuffle  = 1;                      % shuffle random number streams (seed by clock)
cfg.parallel = 1;                      % parallel? true/false
cfg.maxlabs  = 4;                      % maximum number of labs for parallel computing
cfg.nchains  = 4;                      % number of chains
cfg.nswaps = 100000;                    % number of chain swaps
cfg.nsteps = 25;                       % number of steps between chain swaps
cfg.display_status_interval = 5;       % How often to display info
cfg.save_progress_interval = 25;      % How often to save progress 
cfg.adapt_relstep_interval = 100;      % How often to adapt relative step size
cfg.adapt_relstep_rate = 0.24;         % relative step size adaption rate
cfg.optimal_step_acceptance = 0.23;    % optimal rate of step acceptance
cfg.adapt_beta_interval = 250;         % how often to adapt temperature gradient
cfg.adapt_beta_rate = 0.05;            % beta adaption rate
cfg.optimal_swap_acceptance = 0.23;    % optimal rate of swap acceptance
cfg.adapt_last = 20000; %5200;         % last adaption step
cfg.min_adaption_factor = 0.80;        % minimum adaption factor
cfg.max_adaption_factor = 1.25;        % maximum adaption factor
cfg.energy_init_max = 200;             % maximum allowed energy for initialization
cfg.max_beta = 1.0;                    % maximum chain beta (inverse of minimum chain temperature)
cfg.beta_init = 0.667;                 % beta initialization parameter
cfg.relstep_init = 0.075;               % relstep initialization parameter


% setup suffix and regex for progress and init files (DO NOT EDIT)
cfg.progress_suffix = 'progress';
cfg.progress_regex  = sprintf( '%s_%s*.mat', cfg.jobname, cfg.progress_suffix );
cfg.init_suffix = 'init';
cfg.init_regex  = sprintf( '%s_%s*.mat', cfg.jobname, cfg.init_suffix );


%% define custom function handlescfg.data{1}.weight(:,3) = 2;    %Low Dose Akt 308

% initial proposal generator (default is usually ok)
%   template: [params] = @()
cfg.sample_prior_fcn = @() sample_prior(cfg);

% parameter prior (default is usually ok)
%   template: [logpdf] = @(params)
cfg.logpdf_prior_fcn = @(params) logpdf_prior(params,cfg);

% simulation protocols (must provide custom script)
%   template: [err,t,obsv] = @(t,init,params)
cfg.data{1}.protocol_fcn = @(t,init,params,b,c) simulate_pmek(t,init,params,cfg.data{1}.dose,cfg); 
cfg.data{2}.protocol_fcn = @(t,init,params,b,c) simulate_pmek(t,init,params,cfg.data{2}.dose,cfg);
cfg.data{3}.protocol_fcn = @(t,init,params,b,c) simulate_pmek(t,init,params,cfg.data{3}.dose,cfg);
cfg.data{4}.protocol_fcn = @(t,init,params,b,c) simulate_pmek(t,init,params,cfg.data{4}.dose,cfg); 
cfg.data{5}.protocol_fcn = @(t,init,params,b,c) simulate_pmek(t,init,params,cfg.data{5}.dose,cfg);
cfg.data{6}.protocol_fcn = @(t,init,params,b,c) simulate_pmek(t,init,params,cfg.data{6}.dose,cfg);

% universal equilibration protocol (optional)
%   template: [err,t,state] = @(t,init,params)cfg.initial_conditions(1) = 1000;

cfg.equilibrate_fcn = @(t,init,params,b,c) simulate_pmek(t,init,params,cfg.data{1}.dose,cfg);

% heuristic constraint penalty (optional)
%   template: [penalty] = @(d,obsv,params)
% TODO
%cfg.data{1}.heuristic_penalty_fcn = @(d,obsv,params) 0;

% transform simulated data for calibration to experimental data (optional)
%   template: [obsv] = @(obsv,params)
%cfg.transform_sim_for_fit_fcn = @(obsv,params) obsv;

% transform experimental data for plotting graphs (optional)
%   template: [obsv] = @(obsv,params)
%cfg.transform_data_for_plot_fcn = @(obsv,params) exp(obsv);

% energy computation (default is usually ok)
%   template: [energy] = @(params)
cfg.energy_fcn = @(params) energy_generic(params,cfg);

% proposal generator (default is usually ok)
%   template: [params] = @(params,epsilon)
cfg.proposal_fcn = @(params,epsilon) params + epsilon.*randn(1,cfg.nparams);

% update stepsize (default is usually ok)
%   template: [stepsize] = @(relstep)
cfg.update_stepsize_fcn = @(relstep) repmat(relstep, [1 cfg.nparams]).*repmat(cfg.param_scale, [cfg.nchains 1]);

