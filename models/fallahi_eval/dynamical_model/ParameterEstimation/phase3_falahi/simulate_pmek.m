function [err,t,species,obsv] = simulate_params( t, init, params, dose, cfg )

params = 10.^(params);
drug_dose = [3840,12000,37920,120000,379200,1200000,3792000];

if (isempty(t))
    % default timepoints
    t = [cfg.sim_tstart : cfg.sim_dt : cfg.sim_tstop]';
end

if init == 0    %Equil
    
    t = linspace(0,1e5,250)';
    init = zeros(1,53);
    init(1) = params(1);   %grb2
    init(2) = params(2);    %mapk1
    init(3) = params(3);    %kras count
    init(4) = params(4);    %sos1 count
    init(5) = params(5);    %map2k1 count
    init(6) = params(6);    %egfr count
    init(7) = params(7);    %mapk3 count
    init(8) = params(8);    %braf count
    init(9) = params(9);    %egf mRNA count
    init(10) = 0;    %drug  count
    init(11) = params(10);    %phos  count

        params = params(1:61);
       [err, species, obsv] = model_perk(t,init,params);
       
elseif init == 1
    params(1:10) = params(62:71);
    params = params(1:61);
        t = linspace(0,1e5,250)';
    init = zeros(1,53);
    init(1) = params(1);   %grb2
    init(2) = params(2);    %mapk1
    init(3) = params(3);    %kras count
    init(4) = params(4);    %sos1 count
    init(5) = params(5);    %map2k1 count
    init(6) = params(6);    %egfr count
    init(7) = params(7);    %mapk3 count
    init(8) = params(8);    %braf count
    init(9) = params(9);    %egf mRNA count
    init(10) = 0;    %drug  count
    init(11) = params(10);    %phos  count

%         params = params(11:71);
       [err, species, obsv] = model_perk(t,init,params);
    

else
%     dose
	if dose == 1
		init(10) = drug_dose(1);         %Low antigen dose
        % run simulation
                params = params(1:61);
%                 init(10)
        [err, species, obsv] = model_perk(t,init,params);

    elseif dose == 2
		init(10) = drug_dose(4);         %High antigen dose
        % run simulation
                params = params(1:61);
        [err, species, obsv] = model_perk(t,init,params);
        
    elseif dose == 3
		init(10) = drug_dose(7);         %High antigen dose
        % run simulation
                params = params(1:61);
        [err, species, obsv] = model_perk(t,init,params);
	
    end
end

end

