function [params_chain, energy_chain, Percent_Converged] = parallel_tempering( cfg );
%PARALLEL_TEMPERING Generate parameter ensemble using parallel tempering method.
%
%   [params_chain, energy_chain] = parallel_tempering( cfg )
%   generates an ensemble of parameters using the parallel tempering method,
%   where 'cfg' is a configuration struct. See config_clp_ha.m for an example of
%   a script that generates a configuration struct.
%
%   Author: justinshogg@gmail.com
%   Acknowledgments:  S.Lukens, whose code was the starting point for this work.
%   Last Update: 20 Nov 2012
%
%   Features:
%     periodically save progress to file
%     automatic resume if process fails
%     automatic stepsize adjustment
%     automatic temperature adjustment
%     parallel processing
%     stores acceptance rate, temperature and stepsize history
%     modular design (user defined config, log-energy, and proposal generators)
%
%   Assocated files:
%     parallel_temperating.m - main parallel tempering function (this file)
%     config_<model>.m       - configuration file, user defined
%     display_config.m       - display configurations
%     display_chains.m       - display chain status
%     update_stepsize.m      - update stepsize (call if relative_step_size has changed) [this is now obsolete]
%     save_progress.m        - save progress to a file
%     adapt_relstep.m        - adapt relative_step_size to improve step acceptance
%     adapt_beta.m           - adapt beta (temperature) to improve swap acceptance
%     save_init.m            - save an 'init' file that permits chain continuation
%     init_chains.m          - find initial paramters for chains
%
%   CHANGES:
%     13feb2013 - Random number streams are now initialized and shuffled prior to
%                 initializing the chains. Previously, streams were shuffled after
%                 initializing the chains, leading to similar starting points in 
%                 otherwise independent trials. [JSH]
%     17Mar2014 - Adding Gelmin-Rubin Convergence Criteria calculation.
%                 [RPS]
%     26Mar2014 - Adding Bayesian model checking.
%                 [RPS]

nsamples = 0;
naccepted = 0;
fprintf(1,'------------------------------------------------------\n');
fprintf(1,'Parallel Tempering for CLP_HA\n');


% look for progress file
new_chain = 1;
last_swap = 0;
progress_files = dir( cfg.progress_regex );
for file_idx = 1: length(progress_files)
    file = progress_files(file_idx);
    fprintf(1,'------------------------------------------------------\n');
    fprintf(1,'Found progress file %s\n', file.name );
    fprintf(1,'Will attempt to restart aborted chains . . .\n');
    load( file.name );
    new_chain = 0;

    % restarting aborted chain
    % get last params and energt from chain
    params_curr = params_chain(:,:,last_swap+1);
    energy_curr = energy_chain(:,last_swap+1);
    % epsilon = proposal stepsize array
    epsilon = cfg.update_stepsize_fcn( relative_step_size );
    Percent_Converged = converged_params(:,last_swap);
    chisquare = chisquare_value(:,last_swap);

    % what if there is more than one progress file??
    break;
end


% display configuration info
display_config(cfg);


% initialize random number streams
fprintf(1,'------------------------------------------------------\n');
fprintf(1,'Initializing random number streams ');
if (cfg.parallel)
    % start labs
    nlabs = min( [cfg.maxlabs, cfg.nchains] );
    fprintf(1,'[ method=%s, nstreams=%d, shuffle=%d ]\n', 'mrg32k3a', nlabs+1, cfg.shuffle );
    % create independent randstreams
    if (cfg.shuffle)
        [globalstream, labstreams{1:nlabs}] = RandStream.create('mrg32k3a','NumStreams',nlabs+1,'Seed','shuffle');
        %[globalstream, labstreams{1:nlabs}] = RandStream.create('mrg32k3a','NumStreams',nlabs+1,'Seed',1); %version error
    else
        [globalstream, labstreams{1:nlabs}] = RandStream.create('mrg32k3a','NumStreams',nlabs+1,'Seed',1);
    end 

else
    % initialize a single randstream
    fprintf(1,'[ method=%s, nstreams=%d, shuffle=%d ]\n', 'mrg32k3a', 1, cfg.shuffle );
    if (cfg.shuffle)
        globalstream = RandStream('mrg32k3a','Seed','shuffle');
        %globalstream = RandStream('mrg32k3a','Seed',1); %Matlab version error
    else
        globalstream = RandStream('mrg32k3a','Seed',1);
    end
end
% set the global stream now (will set streams for labs later)
RandStream.setGlobalStream( globalstream );
%RandStream.setDefaultStream( globalstream ); %another Matlab version problem

% get proposal generator and energy function
proposal_fcn = cfg.proposal_fcn;
energy_fcn = cfg.energy_fcn;
update_stepsize_fcn = cfg.update_stepsize_fcn;

% start new chain or initialize a restart
if (new_chain)
    % allocate space for working data
    params_curr   = zeros(cfg.nchains, cfg.nparams);       % current params
    energy_curr   = zeros(cfg.nchains, 1);                 % current energy 

    % allocate space for chain data
    params_chain    = zeros(cfg.nchains, cfg.nparams, cfg.nswaps+1);  % parameter samples
    energy_chain    = zeros(cfg.nchains, cfg.nswaps+1);               % energy samples 
    step_acceptance = zeros(cfg.nchains, cfg.nswaps);   % number of accepted steps (per chain per swap)
    swap_acceptance = zeros(cfg.nchains-1, cfg.nswaps); % number of accepted swaps (per chain per swap)
    relstep_history = zeros(cfg.nchains, cfg.nswaps);   % relstep history
    beta_history    = zeros(cfg.nchains, cfg.nswaps);   % beta history
    swap_time       = zeros(cfg.nswaps,1);              % time per swap
    converged_params    = zeros(cfg.nchains, cfg.nswaps);   % percent of parameters converged per G-R
    chisquare_value  = zeros(cfg.nchains, cfg.nswaps);   % chi square goodness of fit

    % look for initialization file
    found_init = 0;
    init_files = dir( cfg.init_regex );
    for file_idx = 1: length(init_files)
        file = init_files(file_idx);
        fprintf(1,'------------------------------------------------------\n');
        fprintf(1,'Found initialization file %s\n', file.name );
        load( file.name );
        found_init = 1;

        % initialize beta and relative step size
        beta                = init.beta;
        relative_step_size  = init.relstep;
        % epsilon = proposal stepsize array
        epsilon = update_stepsize_fcn( relative_step_size );
        % initialize chains
        params_curr = init.params;
        params_chain(:,:,1) = params_curr;
        for chain_idx = 1 : cfg.nchains 
            energy_chain(chain_idx,1) = energy_fcn(params_curr(chain_idx,:,1));
        end
        energy_curr = energy_chain(:,1);
        Percent_Converged = init.converged;
        chisquare_value = init.chisquare;
        % what if there is more than one initialization file??
        clear init;
        break;
    end

    % initialize chains from scratch
    if (~found_init)
        % initialize beta and relative step size
        %beta = ((cfg.beta_init).^[0:(cfg.nchains-1)])';  % beta = inverse chain temperature
        beta = cfg.max_beta * (cfg.beta_init.^[0:(cfg.nchains-1)])';   % beta = inverse chain temperature
        relative_step_size = cfg.relstep_init./beta;      % relative step size w.r.t. log-parameter interval
        % epsilon = proposal stepsize array
        epsilon = update_stepsize_fcn(relative_step_size);
        % initialize chains
        [energy_curr, params_curr] = init_chains( epsilon, cfg );
        params_chain(:,:,1) = params_curr;
        energy_chain(:,1)   = energy_curr;
    end
end


% display additional data
fprintf(1,'------------------------------------------------------\n');
fprintf(1,'Ready to go!\n' );
fprintf(1,'------------------------------------------------------\n');
fprintf(1,'  total swaps = %d\n', cfg.nswaps);
if (last_swap > 0)
    fprintf(1,'  last completed swap = %d\n', last_swap);
end
fprintf(1,'------------------------------------------------------\n');
display_chains( last_swap, cfg, beta, relative_step_size, step_acceptance, swap_acceptance, energy_chain, converged_params, chisquare_value );



% if parallel, start processor pool
if (cfg.parallel)
    fprintf(1,'------------------------------------------------------\n');
    fprintf(1,'Starting Matlab processor pool ');
    % start labs
    nlabs = min( [cfg.maxlabs, cfg.nchains] );
    fprintf(1,'[ nlabs=%d ]\n', nlabs );
    parpool( 'local', nlabs );
    spmd
        % initialize random number stream for each lab
	    RandStream.setGlobalStream( labstreams{labindex} );
        %RandStream.setDefaultStream( globalstream ); %another Matlab version problem
    end
end


% ---------------------------------------------------------------------
%      START CHAINS:
% ---------------------------------------------------------------------
fprintf(1,'------------------------------------------------------\n');
fprintf(1,'Running chains...\n');
counter=1;

% ---- loop over swaps ------------------------------------------------
start_swap = last_swap + 1;
for swap_idx = start_swap : cfg.nswaps

    tic;

    % ---- loop over chains -----------------------------------------------
    n_accept_steps = zeros(cfg.nchains,1);
    parfor chain_idx = 1 : cfg.nchains

        % ---- loop over MCMC steps between swaps -----------------------------
        for step_idx = 1 : cfg.nsteps
    
% %             % get proposed parameters
             params_prop = proposal_fcn( params_curr(chain_idx,:), epsilon(chain_idx,:) );
            
            %Changes by RPS 3/26/14
            %Need energy function to output its 'output' value
            
            % compute proposed energy     
            energy_prop = energy_fcn(params_prop);
            
            % ACCEPT/REJECT
            delta_energy = energy_prop - energy_curr(chain_idx);
            % acceptance probability:
            h = min(1, exp(-beta(chain_idx) * delta_energy));
            if (rand < h)
                % accept proposal!
                params_curr(chain_idx,:) = params_prop;
                energy_curr(chain_idx) = energy_prop;
                step_acceptance(chain_idx,swap_idx) = step_acceptance(chain_idx,swap_idx) + 1; 
            end
                

        % ---------------------------------------------------------------------
        end  % loop over MCMC steps between swaps

    % ---------------------------------------------------------------------
    end  % loop over chains
    

    % Try to Swap Chains (highest temperature to lowest temperature)
    for chain_idx = cfg.nchains : -1 : 2
        
        % Difference.
        delta_energy = energy_curr(chain_idx, 1) - energy_curr(chain_idx-1,1);
        delta_beta = beta(chain_idx) - beta(chain_idx-1);
        
        if ( delta_beta*delta_energy > log(rand(1)) )
            % swap chains . . .
            energy_tmp = energy_curr(chain_idx-1, 1);
            params_tmp = params_curr(chain_idx-1, :);
            
            energy_curr(chain_idx-1, 1) = energy_curr(chain_idx, 1);
            energy_curr(chain_idx, 1)   = energy_tmp;
            
            params_curr(chain_idx-1,:) = params_curr(chain_idx, :);
            params_curr(chain_idx,:) = params_tmp;
            swap_acceptance(chain_idx-1,swap_idx) = 1;
        end
    end

    % ---- check for convergence ------------------------------------------
    if swap_idx < 2000
        Percent_Converged = 0;
    elseif mod(swap_idx,1000) == 0 
%         Percent_Converged = GelmanRubin(swap_idx/2,cfg.nchains,cfg.nparams,params_chain(:,:,1:swap_idx));
        Percent_Converged = 0;

    end
    
    % ---- check for sampling finish ------------------------------------------
    [chisquare] = 0; %chi( params_chain(1,:,swap_idx), cfg );
%     [nsamples,naccepted,outcome] = bmc(nsamples,naccepted,chisquare,cfg);
%     if outcome == 0
% %         break
%     elseif outcome == 1
% %         break
%     elseif outcome == 2
%         nsamples = nsamples + 1;
%     end
            
    
    % ---- save chain -----------------------------------------------------
    params_chain(:,:,swap_idx+1)  = params_curr;
    energy_chain(:,swap_idx+1)    = energy_curr;
    relstep_history(:,swap_idx)   = relative_step_size;
    beta_history(:,swap_idx)      = beta;
    converged_params(:,swap_idx)  = Percent_Converged;
    chisquare_value(:,swap_idx)   = chisquare;

    % remember that this swap was completed
    last_swap = swap_idx;

    % stop timer    
    swap_time(swap_idx) = toc;

    % ---- display output -------------------------------------------------
    if ( mod(swap_idx, cfg.display_status_interval)==0 )
        fprintf(1,'------------------------------------------------------\n');
        fprintf(1,'Finished swap %d\n', swap_idx);
        fprintf(1,'  time for this swap = %.2f sec\n', swap_time(swap_idx) );
        fprintf(1,'  total time = %.2f sec\n', sum(swap_time(1:swap_idx)) );
        display_chains( last_swap, cfg, beta, relative_step_size, step_acceptance, swap_acceptance, energy_chain, converged_params, chisquare_value);
    end

    % ---- adapt step size ------------------------------------------------ 
    if ( mod(swap_idx, cfg.adapt_relstep_interval)==0  &&  mod(swap_idx, cfg.adapt_beta_interval)~=0  &&  swap_idx <= cfg.adapt_last )
        [relative_step_size] = adapt_relstep( relative_step_size, step_acceptance, swap_idx, cfg );
        epsilon = update_stepsize_fcn(relative_step_size);
    end

    % ---- adapt max temperature ----------------------------------------- 
    if ( mod(swap_idx, cfg.adapt_beta_interval)==0  &&  swap_idx <= cfg.adapt_last )
        [beta] = adapt_beta( beta, swap_acceptance, swap_idx, cfg );
    end

    % ---- save progress -------------------------------------------------- 
    if ( mod(swap_idx, cfg.save_progress_interval)==0 )
        save_progress( cfg, swap_idx, params_chain, energy_chain, step_acceptance, swap_acceptance, ...
                       swap_time, last_swap, relative_step_size, beta, relstep_history, beta_history, converged_params, chisquare_value );
    end

%     if Percent_Converged == 1 & mod(swap_idx, cfg.save_progress_interval)==0 %& counter == 1
%         break
%     else 
%         counter = 1;
%     end
end


% ---- display final info ---------------------------------------------
fprintf(1,'------------------------------------------------------\n');
fprintf(1,'Parallel Tempering is complete!\n' );
fprintf(1,'------------------------------------------------------\n');
fprintf(1,'  total swaps = %d\n', last_swap);
fprintf(1,'  total time = %.2f sec\n', sum(swap_time(1:swap_idx)) );
fprintf(1,'------------------------------------------------------\n');
display_chains( last_swap, cfg, beta, relative_step_size, step_acceptance, swap_acceptance, energy_chain, converged_params, chisquare_value );


% ---- save final data ------------------------------------------------
save_progress( cfg, swap_idx, params_chain, energy_chain, step_acceptance, swap_acceptance, ...
               swap_time, last_swap, relative_step_size, beta, relstep_history, beta_history, converged_params, chisquare_value );


% ---- creating init file --------------------------------------------- 
save_init( cfg, params_chain, energy_chain, beta_history, relstep_history, converged_params, chisquare_value );


% ---- close processor pool -------------------------------------------
fprintf(1,'------------------------------------------------------\n');
if (cfg.parallel)
    delete(gcp);
end


% ---- all done ------------------------------------------------------- 




