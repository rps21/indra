function [] = plot_baccam3( params, cfg )
%PLOT_BACCAM3 Plot baccam3 trajectory distribution over parameter ensemble
%
%  [] = plot_baccam3( params, cfg )
%
%  where 'params' is a (N x P) array of parameter samples
%  and 'cfg' is the configuration struct.
 
fig_idx = 1;
fh = figure(fig_idx);
set( fh, 'Color', [1 1 1] );

%% generate trajectories
sim_t = [ cfg.sim_tstart : cfg.sim_dt : cfg.sim_tstop ]';
T = length(sim_t);
O = cfg.nobsv;
N = size(params,1);
sim_obsv_all = zeros(T,O,N);
%energy_all = zeros(N,1);
n_err = 0;

for n = 1:N

    % simulate HA and shammodel_baccam3( (cfg.data{1}.time), y0, params, [],[],[],[], 0 );
    y0 = cfg.initial_conditions;
    y0( cfg.obsv_map.V ) = params(n, cfg.param_map.V0 );
    [err, sim_obsv] = model_baccam3( sim_t, y0, params(n,:), [],[],[],[], 0 );

    % compute energy
    %energy = energy_baccam3( params(n,:), cfg );

    % don't save trajectory if there was an error
    if (~err)
        sim_obsv_all(:,:,n-n_err) = sim_obsv;
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
    sim_obsv_all = sim_obsv_all(:,:,1:(N-n_err));
    %energy_all = energy_all(1:(N-n_err));
end


sim_obsv_all( find(sim_obsv_all < cfg.tolerance) ) = cfg.tolerance;

%% calculate mean and stdev for simulated data
sim_meanlogX   = mean( log10(sim_obsv_all), 3);
sim_stdevlogX  = sqrt( var( log10(sim_obsv_all), [], 3));
%sim_meanX      = 10.^(sim_meanlogX);
%sim_upperbound = 10.^(sim_meanlogX + sim_stdevlogX);
%sim_upperbound( find(isnan(sim_upperbound)) ) = 0;
%sim_lowerbound = 10.^(sim_meanlogX - sim_stdevlogX);
%sim_lowerbound( find(isnan(sim_lowerbound)) ) = 0;
sim_meanX      = (sim_meanlogX);
sim_upperbound = (sim_meanlogX + sim_stdevlogX);
sim_upperbound( find(isnan(sim_upperbound)) ) = 0;
sim_lowerbound = (sim_meanlogX - sim_stdevlogX);
sim_lowerbound( find(isnan(sim_lowerbound)) ) = 0;


% Experimental data
expt_t = cfg.data{1}.time;
expt_meanlogX  = cfg.data{1}.meanlogX;
expt_stdevlogX = cfg.data{1}.stdevlogX;


% plot experimental data and the distribution of model trajectories
for o = 1 : cfg.nobsv

    subplot(1,3,o);
    hold off;
    fill( [sim_t; flipud(sim_t)], [sim_lowerbound(:,o); flipud(sim_upperbound(:,o))], ...
            'b', 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'EdgeAlpha', 0.4);
    hold on;
    plot( sim_t, sim_meanX(:,o), '-r' );
%    errorbar( expt_t, 10.^expt_meanlogX(:,o), ...
%                10.^expt_meanlogX(:,o) - 10.^(expt_meanlogX(:,o) - expt_stdevlogX(:,o)), ...
%                10.^(expt_meanlogX(:,o) + expt_stdevlogX(:,o)) - 10.^expt_meanlogX(:,o), ...
    errorbar( expt_t, expt_meanlogX(:,o),  expt_stdevlogX(:,o),  expt_stdevlogX(:,o), 'sb' );
    xlabel('time', 'fontSize', 12);
    %ylabel( sprintf('log10 %s (%s)', cfg.obsv_names{o}, cfg.obsv_units{o} ), 'fontSize', 12);
    ylabel( sprintf('log %s', cfg.obsv_names{o} ), 'fontSize', 12);
    axis([cfg.sim_tstart cfg.sim_tstop -Inf Inf]);
    hold off;
end
title('baccam3 trajectory distribution','fontsize',14);

