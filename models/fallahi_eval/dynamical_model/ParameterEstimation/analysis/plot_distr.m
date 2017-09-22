function [] = plot_distr( params, cfg, figidx )
%PLOT_DISTR Plot default trajectory distribution over parameter ensemble
%
%  [] = plot_distr( params, cfg )
%  [] = plot_distr( params, cfg, figidx )
%
%  where 'params' is a (S x P) array of parameters and 'cfg' is the
%  configuration struct, and 'figidx' is an optional figure index.
%
%  S = number of samples, P = number of parameters


% permute params if the samples are lie along the 3rd index
if (size(params,1)==1 & size(params,3)>1 )
    params = permute(params,[3 2 1]);
end

% set upper and lower quantiles for distribution fill
upperquant = 0.84;
lowerquant = 0.16;
% number of horizontal subplots
n_horiz_subplots = 4;
% legend location
legend_location = 'NorthEast';


% define xlabel
if isfield(cfg,'time_units')
    xlabel_string = sprintf( 'time (%s)', cfg.time_units );
else
    xlabel_string = 'time';
end
% set up time vector
sim_t = [cfg.sim_tstart : cfg.sim_dt : cfg.sim_tstop]';

% define ylabels
ylabel_strings = {};
for o = 1:cfg.nobsv
    ylabel_strings{o} = sprintf( '%s (%s)', cfg.obsv_names{o}, cfg.obsv_units{o} );
end

% define experiment labels
expt_strings = {};
for d = 1:length(cfg.data)

    if isfield( cfg.data{d}, 'name' )
        expt_strings{d} = cfg.data{d}.name;
    else
        expt_strings{d} = sprintf( 'expr%d', d );
    end
end


% set up linestyles, etc (supports up to three experiments)
linestyle  = {'-b','-g','r','-k'};
fillstyle = {'b','g','r','k'};
errbarstyle = {'o','s','v','x'};
errbarcolor = {[0 0 0.7], [0 0.7 0], [0.7 0 0], [0.1 0.1 0.1]};


% dimensions
T = length(sim_t);     % time points
O = cfg.nobsv;         % observables
D = length(cfg.data);  % experiments
S = size(params,1);    % param samples


% initialize arrays for data storage
for d = 1:D
    sim_obsv{d} = zeros(T,O,S);
end


% loop over parameter samples
for s = 1:S

    % loop over experiments
    for d = 1:D

        % experiment specific settings?
        % equillibrate?

        % simulate w/ this parameter set
        if isfield( cfg.data{d}, 'protocol_fcn' )
            [err, ~, obsv] = cfg.data{d}.protocol_fcn( sim_t, [], params(s,:) );
            if (err) return; end
        else
            [err, ~, obsv] = cfg.data{d}.simulate_fcn( sim_t, [], params(s,:) );
            if (err) return; end
        end

        if isfield(cfg,'transform_data_for_plot_fcn')
            obsv = cfg.transform_data_for_plot_fcn(obsv,params(s,:));
        end

        % save results
        sim_obsv{d}(:,:,s) = obsv;

    end

    % show progress
    fprintf(1,'.');
    if mod(s, 72)==0
        fprintf(1,'\n');
    end

end
fprintf(1,'\n');


%% calculate mean and stdev for simulated data
for d = 1:D
    sim_mean{d} = quantile(sim_obsv{d}, 0.5, 3);
    sim_upperbound{d} = quantile(sim_obsv{d}, upperquant, 3);
    sim_lowerbound{d} = quantile(sim_obsv{d}, lowerquant, 3);
end


% fetch experimental data
expt = cfg.data;

% plot experimental data and the distribution of model trajectories
if exist('figidx')  fh = figure(figidx);
else  fh = figure;  end
set( fh, 'Color', [1 1 1] );
for o = 1:cfg.nobsv
    
    subplot( ceil(cfg.nobsv/n_horiz_subplots), min(n_horiz_subplots,cfg.nobsv), o);
    hold off;
    for d = D:-1:1

        fill( [sim_t; flipud(sim_t)], [sim_lowerbound{d}(:,o); flipud(sim_upperbound{d}(:,o))], ...
                fillstyle{d}, 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'EdgeAlpha', 0.4);
        hold on;
        plot( sim_t, sim_mean{d}(:,o), linestyle{d} );

        % plot errorbars (if any)
        if ( any(~isnan(expt{d}.mean(:,o))) )
            ebh = errorbar( expt{d}.time , expt{d}.mean(:,o), ...
                        expt{d}.stdev(:,o), expt{d}.stdev(:,o), ...
                        errbarstyle{d}, 'color', errbarcolor{d}, 'linewidth', 1.1 );
            % increase width of errorbars
            adjust_errorbars(ebh, 0.2);
        end
       
    end

    if (o==cfg.nobsv)
        lh = legend( expt_strings );
        set(lh,'Location',legend_location);
    end

    % set axes limits
    if ( isfield(cfg,'plot_xlim') )
        xlim( cfg.plot_xlim );
    else
        xlim( [cfg.sim_tstart cfg.sim_tstop] );
    end
    if ( isfield(cfg,'plot_ylims') )
        ylim( cfg.plot_ylims{o} );
    else
        ylim( [-Inf Inf] );
    end
    % label axes
    xlabel( xlabel_string, 'fontSize', 12);
    ylabel( ylabel_strings{o}, 'fontSize', 12);

    hold off;
end

% end plot_distr function
end



function [] = adjust_errorbars( ebh, scale )
    % increase width of errorbars
    ebh = get(ebh,'children');
    bardata = get(ebh(2),'Xdata');
    left_idx = sort( [4:9:length(bardata), 7:9:length(bardata)] );
    right_idx = sort( [5:9:length(bardata), 8:9:length(bardata)] );
    bardata(left_idx)  = bardata(left_idx)  - scale;
    bardata(right_idx) = bardata(right_idx) + scale;
    set(ebh(2),'Xdata',bardata);
end

