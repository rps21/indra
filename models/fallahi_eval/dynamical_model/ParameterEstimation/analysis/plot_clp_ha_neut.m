function [] = plot_clp_ha( params, cfg  )
%PLOT_CLP_HA Plot CLP trajectory distribution over parameter ensemble for HA and Sham
%
%  [] = plot_clp_ha( params, cfg, fig_idx )
%
%  where 'params' is a (N x P) array of parameters and 'cfg' is the configuration struct.
%  N = number of samples, P = number of parameters
 

fh = figure; %(fig_idx);
set( fh, 'Color', [1 1 1] );

%% generate trajectories
sim_t_sh = [ cfg.sim_tstart : cfg.sim_dt : cfg.sim_tstop ]';
sim_t_ha = [ cfg.sim_tstart : cfg.sim_dt : cfg.sim_tstop ]';
T_sh = length(sim_t_sh);
T_ha = length(sim_t_ha);
D = cfg.nobsv;
N = size(params,1);
sim_obsv_sh_all = zeros(T_sh,D,N);
sim_obsv_ha_all = zeros(T_ha,D,N);
%energy_all = zeros(N,1);
n_err = 0;

for n = 1:N

    % simulate HA and sham
    [err_sh, sim_obsv_sh] = simulate_neutrophil( sim_t_sh, [], params(n,:), 0, cfg );
    [err_ha, sim_obsv_ha] = simulate_neutrophil( sim_t_ha, [], params(n,:), 1, cfg );

    % compute energy
    %energy = energy_clp_ha( params(n,:), cfg );

    % don't save trajectory if there was an error
    if (~err_sh && ~err_ha)
        sim_obsv_sh_all(:,:,n-n_err) ...
            = correct_for_dilution( counts_to_conc( sim_obsv_sh, params(n,:), cfg), params(n,:), cfg);
        sim_obsv_ha_all(:,:,n-n_err) ...
            = correct_for_dilution( counts_to_conc( sim_obsv_ha, params(n,:), cfg), params(n,:), cfg);;
        %energy_all(n-n_err) = energy;
    else
        n_err = n_err+1;
    end

    if mod(n,5)==0
        fprintf(1,'.');
        if mod(n, 5*79)==0
            fprintf(1,'\n');
        end
    end
end
fprintf(1,'\n');

% trim empty trajectoires
if ( n_err )
    sim_obsv_sh_all = sim_obsv_sh_all(:,:,1:(N-n_err));
    sim_obsv_ha_all = sim_obsv_ha_all(:,:,1:(N-n_err));
    %energy_all = energy_all(1:(N-n_err));
end



%% calculate mean and stdev for simulated data
% Sham
sim_t_sh = sim_t_sh/cfg.mins_per_hr;
sim_sh_meanlogX   = mean( log2(sim_obsv_sh_all), 3);
sim_sh_stdevlogX  = sqrt( var( log2(sim_obsv_sh_all), [], 3));
sim_sh_meanX      = 2.^(sim_sh_meanlogX);
sim_sh_upperbound = 2.^(sim_sh_meanlogX + sim_sh_stdevlogX);
sim_sh_upperbound( find(isnan(sim_sh_upperbound)) ) = 0;
sim_sh_lowerbound = 2.^(sim_sh_meanlogX - sim_sh_stdevlogX);
sim_sh_lowerbound( find(isnan(sim_sh_lowerbound)) ) = 0;
% HA
sim_t_ha = sim_t_ha/cfg.mins_per_hr;
sim_ha_meanlogX  = mean( log2(sim_obsv_ha_all), 3);
sim_ha_stdevlogX = sqrt( var( log2(sim_obsv_ha_all), [], 3));
sim_ha_meanX      = 2.^(sim_ha_meanlogX);
sim_ha_upperbound = 2.^( sim_ha_meanlogX + sim_ha_stdevlogX );
sim_ha_upperbound( find(isnan(sim_ha_upperbound)) ) = 0;
sim_ha_lowerbound = 2.^( sim_ha_meanlogX - sim_ha_stdevlogX );
sim_ha_lowerbound( find(isnan(sim_ha_lowerbound)) ) = 0;


% Sham
expt_t_sh = cfg.data{1}.time;
expt_sh_meanlogX  = cfg.data{1}.meanlogX;
expt_sh_stdevlogX = cfg.data{1}.stdevlogX;
% HA
expt_t_ha = cfg.data{2}.time;
expt_ha_meanlogX  = cfg.data{2}.meanlogX;
expt_ha_stdevlogX = cfg.data{2}.stdevlogX;

% plot experimental data and the distribution of model trajectories
for s = 1 : cfg.nobsv

    subplot(2,4,s);
    hold off;
    fill( [sim_t_sh; flipud(sim_t_sh)], [sim_sh_lowerbound(:,s); flipud(sim_sh_upperbound(:,s))], ...
            'r', 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'EdgeAlpha', 0.4);
    hold on;
    fill( [sim_t_ha; flipud(sim_t_ha)], [sim_ha_lowerbound(:,s); flipud(sim_ha_upperbound(:,s))], ...
            'b', 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'EdgeAlpha', 0.4);

    plot( sim_t_sh, sim_sh_meanX(:,s), '-r' );
    plot( sim_t_ha, sim_ha_meanX(:,s), '-b' );
    errorbar( expt_t_sh, 2.^expt_sh_meanlogX(:,s), ...
                2.^expt_sh_meanlogX(:,s) - 2.^(expt_sh_meanlogX(:,s) - expt_sh_stdevlogX(:,s)), ...
                2.^(expt_sh_meanlogX(:,s) + expt_sh_stdevlogX(:,s)) - 2.^expt_sh_meanlogX(:,s), ...
                'sr' );
    errorbar( expt_t_ha, 2.^expt_ha_meanlogX(:,s), ...
                2.^expt_ha_meanlogX(:,s) - 2.^(expt_ha_meanlogX(:,s) - expt_ha_stdevlogX(:,s)), ...
                2.^(expt_ha_meanlogX(:,s) + expt_ha_stdevlogX(:,s)) - 2.^expt_ha_meanlogX(:,s), ...
                'sb' );
    xlabel('time (hrs)', 'fontSize', 12);
    ylabel( sprintf('%s (%s)', cfg.obsv_names{s}, cfg.obsv_units{s} ), 'fontSize', 12);
    axis([0 cfg.sim_tstop/cfg.mins_per_hr 0 Inf]);
    hold off;
end

