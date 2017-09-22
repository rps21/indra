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
N = size(params,2);
%N = size(params,1);
sim_obsv_all = zeros(T,O,N);
sim_obsv_all2 = zeros(T,O,N);

    params = 10.^(params);


%energy_all = zeros(N,1);
n_err = 0;

for n = 1:N
%     % simulate HA and shammodel_baccam3( (cfg.data{1}.time), y0, params, [],[],[],[], 0 );
    y0 = cfg.initial_conditions; %low dose
        y0(1) = params(6);
        y0(2) = params(7);
        y0(3) = params(8);
        y0(4) = params(9);
        y0(5) = params(10);
        y0(6) = params(11);
        y0(7) = params(12);
%     [err, sim_obsv] = model_IL17( sim_t, y0, params(:,n)', [],[],[],[],0 ); %(t,init,params,[],[],[],[],0)
    [err, sp, sim_obsv] = model_IL17(sim_t, y0, params(:,n)');



    
    % compute energy
    energy = energy_generic( params(:,n)', cfg );

    % don't save trajectory if there was an error
    if (~err)
        sim_obsv_all(:,:,n-n_err) = sim_obsv;

%         energy_all(n-n_err) = energy;
    else
        n_err = n_err+1
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
sim_obsv_all2( find(sim_obsv_all2 < cfg.tolerance) ) = cfg.tolerance;

% % calculate mean and stdev for simulated data
% sim_meanlogX   = mean( log10(sim_obsv_all), 3);
% sim_stdevlogX  = sqrt( var( log10(sim_obsv_all), [], 3));
% sim_meanX      = 10.^(sim_meanlogX);
% sim_upperbound = 10.^(sim_meanlogX + sim_stdevlogX);
% sim_upperbound( find(isnan(sim_upperbound)) ) = 0;
% sim_lowerbound = 10.^(sim_meanlogX - sim_stdevlogX);
% sim_lowerbound( find(isnan(sim_lowerbound)) ) = 0;
% sim_meanX      = (sim_meanlogX);
% sim_upperbound = (sim_meanlogX + sim_stdevlogX);
% sim_upperbound( find(isnan(sim_upperbound)) ) = 0;
% sim_lowerbound = (sim_meanlogX - sim_stdevlogX);
% sim_lowerbound( find(isnan(sim_lowerbound)) ) = 0;

%% calculate mean and stdev for simulated data - now with no logs
sim_meanlogX   = mean( sim_obsv_all, 3);
sim_stdevlogX  = sqrt( var( sim_obsv_all, [], 3));

% sim_meanlogX   = mean( sim_obsv_all);
% sim_stdevlogX  = sqrt( var( sim_obsv_all, []));
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



% % Experimental data
expt_t = cfg.data{1}.time(1:9);
% expt_t = cfg.data{1}.time(1:8);
expt_mean  = cfg.data{1}.mean;
expt_stdev = cfg.data{1}.stdev;

% plot experimental data and the distribution of model trajectories
for o = 1:2 %cfg.nobsv
% 
    subplot(1,2,o);
%     hold off;
% % 
%     fill( [sim_t; flipud(sim_t)], [sim_lowerbound(:,o); flipud(sim_upperbound(:,o))], ...
%            'b', 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'EdgeAlpha', 0.4);
    hold on;
%     plot( sim_t, sim_meanX(:,o), '-r','LineWidth',2 );

    for i = 1:5
        plot(sim_t, sim_obsv_all(:,o,i), '-r')
    end

%    errorbar( expt_t, 10.^expt_meanlogX(:,o), ...
%                10.^expt_meanlogX(:,o) - 10.^(expt_meanlogX(:,o) - expt_stdevlogX(:,o)), ...
%                10.^(expt_meanlogX(:,o) + expt_stdevlogX(:,o)) - 10.^expt_meanlogX(:,o), ...

    errorbar( expt_t, expt_mean(:,o),  expt_stdev(:,o),  expt_stdev(:,o), 'sb' );
    xlabel('Time (seconds)', 'fontSize', 14);
    %ylabel( sprintf('log10 %s (%s)', cfg.obsv_names{o}, cfg.obsv_units{o} ), 'fontSize', 12);
    ylabel( sprintf( cfg.obsv_names{o} ), 'fontSize', 14); %
    axis([cfg.sim_tstart cfg.sim_tstop -Inf Inf]);
     hFig = figure(1);
    set(gcf,'PaperPositionMode','auto')
    set(hFig, 'Position', [0 0 851 293])
    hold off;
end
%title('baccam3 trajectory distribution','fontsize',14);

