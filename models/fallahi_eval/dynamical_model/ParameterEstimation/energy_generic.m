function [energy] = energy_generic( params, cfg )
%ENERGY_GENERIC Calculate parameter energy (log-likelihood) for generic model
% 
%  [energy] = energy_generic( params, cfg )

% check prior probability on parameters
logprior = cfg.logpdf_prior_fcn(params);
if isinf(logprior)
    energy = cfg.big_energy;
    return;
end
energy = -logprior;

% start integration timer
simtimer = tic;

% equillibrate, if required
if isfield(cfg, 'equilibrate_fcn')
    [err,~,sp,obsv_equil] = cfg.equilibrate_fcn( [], 0, params, 1, cfg );
    init1 = sp(end,:);
     
     [err,~,sp,obsv_equil] = cfg.equilibrate_fcn( [], 1, params, 1, cfg );
    init2 = sp(end,:);   
    
    if (err)
        energy = cfg.big_energy;

        return;
    end
else
    init = [];
end

% simulate experiments   
for d = 6%1:2%cfg.nexpt

    if d <= 3
        init = init1;
    else
        init = init2;
    end
    
    % simulate experiment (default initial conditions)
    [err,~,~,obsv_raw] = cfg.data{d}.protocol_fcn( [], init, params, d, cfg );

    %11_12 - covert obsv to smaller array using time map
    fields = fieldnames(cfg.data{1}.time_map);

    for obsvcount = 1:cfg.nobsv
        for timepoint = 1:numel(fields)
%             cfg.data{1}.time_map.(fields{timepoint})
            obsv(timepoint,obsvcount) = obsv_raw(cfg.data{1}.time_map.(fields{timepoint}),obsvcount); %Gives same timepoints as data
        end
    end

    %Model observables
% cfg.obsv_defs = { ...
% 1  struct('name','p-mek', 'units','log10'), ...
% 2  struct('name','p-erk', 'units','log10'), ...
% 3  struct('name','akt_p308', 'units','log10'), ...
% 4  struct('name','akt_p473', 'units','log10'), ...
% 5  struct('name','p-mtor', 'units','log10'), ...
% 6  struct('name','p-s6k', 'units','log10'), ...
% 7  struct('name','p-s6', 'units','log10'), ...
% 8  struct('name','p-p27', 'units','log10'), ...
    
    if obsv(1,1) < .0011
        obsv(1,1) = .0011;
    end
        
    obsv(:,1) = log2(obsv(:,1)./obsv(1,1));   
    obsv = obsv(:,1);
%     obsv(:,2) = log2(obsv(:,2)./obsv(1,2));         
%     obsv(:,3) = log2(obsv(:,3)./obsv(1,3));         
%     obsv(:,4) = log2(obsv(:,4)./obsv(1,4));         
%     obsv(:,5) = log2(obsv(:,5)./obsv(1,5));         
%     obsv(:,6) = log2(obsv(:,6)./obsv(1,6));         
%     obsv(:,7) = log2(obsv(:,7)./obsv(1,7));         
%     obsv(:,8) = log2(obsv(:,8)./obsv(1,8));         

    
    
%     obsv(:,1) = obsv(:,1)./obsv(1,1);         %A20 mRNA
%     obsv(:,2) = obsv(:,2)./obsv(1,2);         %A20 protein
%     obsv(:,3) = obsv(:,3)./obsv(1,3);         %IkB mRNA
%     obsv(:,4) = obsv(:,4)./obsv(1,4);         %IkB protein
%     obsv(:,5) = obsv(:,5)./max(obsv(:,5));    %TAK1~p
     
    if (err)
        energy = cfg.big_energy;
        return;
    end

    % heuristic penalties
    if isfield(cfg.data{d}, 'heuristic_penalty_fcn')
        penalty = cfg.data{d}.heuristic_penalty_fcn(obsv, params);
        if isinf(penalty)
            energy = cfg.big_energy;
            return;
        end
        energy = energy + penalty;
    end

    % if necessary, transform simulated trajectory for computing fitness
    if isfield(cfg, 'transform_sim_for_fit_fcn')
        obsv = cfg.transform_sim_for_fit_fcn(obsv,param);
    end

    % calculate log-likelihood as weighted sum of square errors  
%     obsv
% (obsv - cfg.data{d}.mean).^2
    loglike = nansum(nansum( -cfg.data{1}.weight .* (obsv - cfg.data{d}.mean).^2 ./ (2*(cfg.data{d}.stdev).^2) )); %%%Changed 1 to d, using different mean for each experiment, as each different initial condition
    
% subtract likelihood from energy
    energy = energy - loglike;

end

% penalize for slow integrations
dt = toc(simtimer);
energy = energy + cfg.timepenalty*dt^2;

% all done
return;
