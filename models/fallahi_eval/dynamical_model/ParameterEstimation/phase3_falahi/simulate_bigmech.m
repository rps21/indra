function [err,t,species,obsv] = simulate_params( t, init, params, dose, cfg )

params = 10.^(params);

if (isempty(t))
    % default timepoints
    t = [cfg.sim_tstart : cfg.sim_dt : cfg.sim_tstop]';
end

if init == 0    %Equil
    
    t = linspace(0,1e5,250)';
    init = zeros(1,725);
    init(1) = params(72);
    init(2) = params(55);    %rptor
    init(3) = params(56);    %raf1 count
    init(4) = params(57);    %grb2 count
    init(5) = params(58);    %mapk1 count
    init(6) = params(59);    %sos1 count
    init(7) = params(60);    %mtor count
    init(8) = params(61);    %rps6kb1 count
    init(9) = params(62);    %mapk3 mRNA count
    init(10) = params(63);    %braf  count
    init(11) = params(64);    %akt count
    init(12) = params(65);    %egfr mRNA count
    init(13) = params(66);    %kras count
    init(14) = params(67);    %rictor count
    init(15) = params(68);    %map2k1 count
    init(16) = params(69);    %cdkn1b count    
        init(17) = params(70);    %s6 count    
    init(18) = 0;    %drug count    

       [err, species, obsv] = model_bigmech_test(t,init,params);

else
	if dose == 1
		init(18) = params(71);         %Low antigen dose
%         init(2) = params(2);
        % run simulation
        [err, species, obsv] = model_bigmech_test(t,init,params);

%     elseif dose == 2
% 		init(1) = params(1)*4;         %High antigen dose
%         init(2) = params(2);
%         % run simulation
%         [err, species, obsv] = model_tcell_v26(t,init,params);
	
    end
end

end

