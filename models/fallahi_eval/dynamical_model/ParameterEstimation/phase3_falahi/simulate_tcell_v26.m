function [err,t,species,obsv] = simulate_tcell_v25( t, init, params, dose, cfg )

params = 10.^(params);

if (isempty(t))
    % default timepoints
    t = [cfg.sim_tstart : cfg.sim_dt : cfg.sim_tstop]';
end

if init == 0    %Equil
    
    t = linspace(0,1e8,500)';
    init = zeros(1,66);

    init(1) = 0;
    init(2) = 0;    %Anti CD28
    init(3) = params(3);    %TCR count
    init(4) = params(3);    %CD28 count
    init(5) = params(4);    %PIP count
    init(6) = params(5);    %AKT count
    init(7) = params(6);    %FOXO1 count
    init(8) = params(7);    %PTEN count
    init(9) = params(8);    %PTEN mRNA count
    init(10) = params(9);    %mTORC2  count
    init(11) = params(10);    %NEDD4 count
    init(12) = params(11);    %CK2 mRNA count
    init(13) = params(12);    %PDK1 count
    init(14) = params(13);    %PI3K count
    init(15) = params(14);    %PP2A count
    init(16) = params(15);    %DNA count

    [err, species, obsv] = model_tcell_v26(t,init,params);

else
	if dose == 1
		init(1) = params(1);         %Low antigen dose
        init(2) = params(2);
        % run simulation
        [err, species, obsv] = model_tcell_v26(t,init,params);

    elseif dose == 2
		init(1) = params(1)*4;         %High antigen dose
        init(2) = params(2);
        % run simulation
        [err, species, obsv] = model_tcell_v26(t,init,params);
	
    end
end

end

