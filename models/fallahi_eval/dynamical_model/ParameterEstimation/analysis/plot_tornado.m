function [] = plot_tornado( params, p1, p2, cfg );
%PLOT_TORNADO Plot pairs of parmeters from multiple chains with samples along the z-axis
%
%  [] = plot_tornado( params, p1, p2, cfg )
%
%  where 'params' is a (C x P x N) array of parameter samples,
%  with C=number of chains, P=number of parametersz, N=number of samples;
%  'p1' is the index of the first parameter, 'p2' is the index of the second
%  and 'cfg' is the configuration struct


nchains = size(params,1);
nparams = size(params,2);
nsamples = size(params,3);

style = { 'or', 'oy', 'og', 'ob', 'oc', 'dr', 'dy', 'dg', 'db', 'dc' };

fh = figure; %(1);
set( fh, 'Color', [1 1 1] );

T = [1:nsamples]';
for c = nchains:-1:1
    scatter3( permute(params(c,p1,:), [3 2 1]), permute(params(c,p2,:), [3 2 1]), T, 9, style{c} );
    hold on;
end
hold off;

% axis( [cfg.param_range(1,p1), cfg.param_range(2,p1), cfg.param_range(1,p2), cfg.param_range(2,p2), 0, nsamples] );
xlabel( sprintf('%s', cfg.param_names{p1}) );
ylabel( sprintf('%s', cfg.param_names{p2}) );
zlabel( 'samples' );

