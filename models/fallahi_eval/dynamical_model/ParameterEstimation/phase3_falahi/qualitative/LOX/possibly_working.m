function [err, timepoints, species_out, observables_out ] = possibly_working( timepoints, species_init, parameters, suppress_plot )
%POSSIBLY_WORKING Integrate reaction network and plot observables.
%   Integrates the reaction network corresponding to the BioNetGen model
%   'possibly_working' and then (optionally) plots the observable trajectories,
%   or species trajectories if no observables are defined. Trajectories are
%   generated using either default or user-defined parameters and initial
%   species values. Integration is performed by the CVode library interfaced
%   to MATLAB via the MEX interface. Before running this script, the model
%   source in file possibly_working_cvode.c must be compiled (see that file for details).
%   POSSIBLY_WORKING returns an error value, a vector of timepoints,
%   species trajectories, and observable trajectories.
%   
%   [err, timepoints, species_out, observables_out]
%        = possibly_working( timepoints, species_init, parameters, suppress_plot )
%
%   INPUTS:
%   -------
%   timepoints      : column vector of time points returned by integrator.
%   species_init    : row vector of 53 initial species populations.
%   parameters      : row vector of 62 model parameters.
%   suppress_plot   : 0 if a plot is desired (default), 1 if plot is suppressed.
%
%   Note: to specify default value for an input argument, pass the empty array.
%
%   OUTPUTS:
%   --------
%   err             : 0 if the integrator exits without error, non-zero otherwise.
%   timepoints      : a row vector of timepoints returned by the integrator.
%   species_out     : array of species population trajectories
%                        (columns correspond to species, rows correspond to time).
%   observables_out : array of observable trajectories
%                        (columns correspond to observables, rows correspond to time).
%
%   QUESTIONS about the BNG Mfile generator?  Email justinshogg@gmail.com



%% Process input arguments

% define any missing arguments
if ( nargin < 1 )
    timepoints = [];
end

if ( nargin < 2 )
    species_init = [];
end

if ( nargin < 3 )
    parameters = [];
end

if ( nargin < 4 )
    suppress_plot = 0;
end


% initialize outputs (to avoid error msgs if script terminates early
err = 0;
species_out     = [];
observables_out = [];


% setup default parameters, if necessary
if ( isempty(parameters) )
   parameters = [ 1, 5.000000e+04, 5.000000e+04, 5.000000e+04, 5.000000e+04, 75000, 1.000000e+04, 5.000000e+04, 150000, 3750, 1e5, 1.000000e-06, 1.000000e-02, 1.000000e-06, 1.000000e-01, 1.000000e-06, 1.000000e-01, 1.000000e-06, 1.000000e-01, 1.000000e-06, 1.000000e-06, 1.000000e-01, 1.000000e+02, 1.000000e-06, 1.000000e-01, 1.000000e+02, 1.000000e-06, 1.000000e-01, 1.000000e+02, 1.000000e-06, 1.000000e-01, 1.000000e+02, 1.000000e-06, 1.000000e-01, 1.000000e+02, 1.000000e-06, 1.000000e-01, 1.000000e+02, 1.000000e-06, 1.000000e-01, 1.000000e+02, 1.000000e-06, 1.000000e-01, 1.000000e+02, 1.000000e-06, 1.000000e-01, 1.000000e+02, 1.000000e-06, 1.000000e-01, 1.000000e+02, 1.000000e-06, 1.000000e-01, 1.000000e+02, 1e-6, 1e-1, 1e-5, 5e-6, 1e-1, 10, 5e-7, 1e-1, 5e-4 ];
end
% check that parameters has proper dimensions
if (  size(parameters,1) ~= 1  ||  size(parameters,2) ~= 62  )
    fprintf( 1, 'Error: size of parameter argument is invalid! Correct size = [1 62].\n' );
    err = 1;
    return;
end

% setup default initial values, if necessary
if ( isempty(species_init) )
   species_init = initialize_species( parameters );
end
% check that species_init has proper dimensions
if (  size(species_init,1) ~= 1  ||  size(species_init,2) ~= 53  )
    fprintf( 1, 'Error: size of species_init argument is invalid! Correct size = [1 53].\n' );
    err = 1;
    return;
end

% setup default timepoints, if necessary
if ( isempty(timepoints) )
   timepoints = linspace(0,10,20+1)';
end
% check that timepoints has proper dimensions
if (  size(timepoints,1) < 2  ||  size(timepoints,2) ~= 1  )
    fprintf( 1, 'Error: size of timepoints argument is invalid! Correct size = [t 1], t>1.\n' );
    err = 1;
    return;
end

% setup default suppress_plot, if necessary
if ( isempty(suppress_plot) )
   suppress_plot = 0;
end
% check that suppress_plot has proper dimensions
if ( size(suppress_plot,1) ~= 1  ||  size(suppress_plot,2) ~= 1 )
    fprintf( 1, 'Error: suppress_plots argument should be a scalar!\n' );
    err = 1;
    return;
end

% define parameter labels (this is for the user's reference!)
param_labels = { 'drug_0', 'GRB2_0', 'MAPK1_0', 'KRAS_0', 'SOS1_0', 'MAP2K1_0', 'EGFR_0', 'MAPK3_0', 'BRAF_0', 'EGF_0', 'phos_0', 'kf_kb_act_1', 'kf_e_autophos_1', 'kf_eg_bind_1', 'kr_eg_bind_1', 'kf_ee_bind_1', 'kr_ee_bind_1', 'kf_gs_bind_1', 'kr_gs_bind_1', 'kf_sk_gef_1', 'kf_bm_bind_1', 'kr_bm_bind_1', 'kc_bm_phosphorylation_1', 'kf_bm_bind_2', 'kr_bm_bind_2', 'kc_bm_phosphorylation_2', 'kf_bm_bind_3', 'kr_bm_bind_3', 'kc_bm_phosphorylation_3', 'kf_mm_bind_1', 'kr_mm_bind_1', 'kc_mm_phosphorylation_1', 'kf_mm_bind_2', 'kr_mm_bind_2', 'kc_mm_phosphorylation_2', 'kf_mb_bind_1', 'kr_mb_bind_1', 'kc_mb_phosphorylation_1', 'kf_mm_bind_3', 'kr_mm_bind_3', 'kc_mm_phosphorylation_3', 'kf_mm_bind_4', 'kr_mm_bind_4', 'kc_mm_phosphorylation_4', 'kf_mb_bind_2', 'kr_mb_bind_2', 'kc_mb_phosphorylation_2', 'kf_mm_bind_5', 'kr_mm_bind_5', 'kc_mm_phosphorylation_5', 'kf_mm_bind_6', 'kr_mm_bind_6', 'kc_mm_phosphorylation_6', 'kf_drug_braf', 'kr_drug_braf', 'mapk1_p_decay_slow', 'kf_erk_phos', 'kr_erk_phos', 'kc_erk_phos', 'kf_drug_phos', 'kr_drug_phos', 'map2k1_decay' };



%% Integrate Network Model
try 
    % run simulation
    [err, species_out, observables_out] = possibly_working_cvode( timepoints, species_init, parameters );
catch
    fprintf( 1, 'Error: some problem integrating ODE network! (CVODE exitflag %d)\n', err );
    err = 1;
    return;
end



%% Plot Output, if desired

if ( ~suppress_plot )
    
    % define plot labels
    observable_labels = { 'erk_phos', 'totalerk', 'active_map2k1', 'mek_erk', 'totalbraf', 'active_braf', 'phos_free', 'phos_drugged', 'drug_total', 'drug_bound' };

    % construct figure
    plot(timepoints,observables_out);
    title('possibly_working observables','fontSize',14,'Interpreter','none');
    axis([0 timepoints(end) 0 inf]);
    legend(observable_labels,'fontSize',10,'Interpreter','none');
    xlabel('time','fontSize',12,'Interpreter','none');
    ylabel('number or concentration','fontSize',12,'Interpreter','none');

end



%~~~~~~~~~~~~~~~~~~~~~%
% END of main script! %
%~~~~~~~~~~~~~~~~~~~~~%



% initialize species function
function [species_init] = initialize_species( params )

    species_init = zeros(1,53);
    species_init(1) = params(2);
    species_init(2) = params(3);
    species_init(3) = params(4);
    species_init(4) = params(5);
    species_init(5) = params(6);
    species_init(6) = params(7);
    species_init(7) = params(8);
    species_init(8) = params(9);
    species_init(9) = params(10);
    species_init(10) = 0;
    species_init(11) = params(11);
    species_init(12) = 0;
    species_init(13) = 0;
    species_init(14) = 0;
    species_init(15) = 0;
    species_init(16) = 0;
    species_init(17) = 0;
    species_init(18) = 0;
    species_init(19) = 0;
    species_init(20) = 0;
    species_init(21) = 0;
    species_init(22) = 0;
    species_init(23) = 0;
    species_init(24) = 0;
    species_init(25) = 0;
    species_init(26) = 0;
    species_init(27) = 0;
    species_init(28) = 0;
    species_init(29) = 0;
    species_init(30) = 0;
    species_init(31) = 0;
    species_init(32) = 0;
    species_init(33) = 0;
    species_init(34) = 0;
    species_init(35) = 0;
    species_init(36) = 0;
    species_init(37) = 0;
    species_init(38) = 0;
    species_init(39) = 0;
    species_init(40) = 0;
    species_init(41) = 0;
    species_init(42) = 0;
    species_init(43) = 0;
    species_init(44) = 0;
    species_init(45) = 0;
    species_init(46) = 0;
    species_init(47) = 0;
    species_init(48) = 0;
    species_init(49) = 0;
    species_init(50) = 0;
    species_init(51) = 0;
    species_init(52) = 0;
    species_init(53) = 0;

end


end
