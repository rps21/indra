function [cfg] = init_observable_defs( obsv_defs, cfg )
%INIT_OBSERVABLE_DEFS read observable definitions and set up structs
    cfg.nobsv = length(obsv_defs);
    cfg.obsv_names = {};
    cfg.obsv_units = {};
    for o = 1:cfg.nobsv
        % map from observable names to indices
        cfg.obsv_map.(obsv_defs{o}.name) = o;
        % set up structs
        cfg.obsv_names{o} = obsv_defs{o}.name;
        cfg.obsv_units{o} = obsv_defs{o}.units;
    end

end
