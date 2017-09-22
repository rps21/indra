function [] = plot_hist_energy( samples, energy, cfg );
%PLOT_HIST_ENERGY plot histograms of parameter samples w/ energy
%
%  [] = plot_histogram( samples, energy, cfg )
%  
%  where 'samples' is a (C x P x N) array of parameter samples,
%  'energy' is a (C x N) array of energy value for the samples
%  with C=number of chains, P=number of parameters, N=number of samples;
%  and 'cfg' is the configuration struct.
%
%  Histograms are shown for one parameter at a time. Hit any key to see
%  the next histogram.

%added first to lines from plot_histogram
width = 1.6;
param_range = [ cfg.param_location - width*cfg.param_scale; cfg.param_location + width*cfg.param_scale ];
%param_range = cfg.param_range;
samples = permute(samples, [3 2 1]);

nsamples = size(samples,1);
nparams = size(samples,2);
nchains = size(samples,3);

fh = figure;
set( fh, 'Color', [1 1 1] );

res=40;
for p = 1:nparams

    if param_range(1,p) == param_range(2,p)
        continue
    end
    
    edges = linspace( param_range(1,p), param_range(2,p), res )';

    maxpdf = 0;
    for c = 1:nchains

        [N] = histc( samples(:,p,c), edges(1:end-1) );
        pdf = N/sum((edges(2:end)-edges(1:end-1)).*N);
        maxpdf = max( [maxpdf; pdf ] );
        subplot(2,nchains,c);
        bar(edges(1:end-1),pdf,'histc');
        xlabel( sprintf('%s', cfg.param_names{p}) );
        ylabel( 'pdf' );
        title( sprintf('chain %d', c) );

        bin_center = zeros(length(pdf)-1,1);
        avg_energy = zeros(length(pdf)-1,1);
        for k = 1 : length(pdf)-1
            bink = find( and(edges(k) <= samples(:,p,c), samples(:,p,c) < edges(k+1)) );
            avg_energy(k) = mean(energy(c,bink));
            bin_center(k) = (edges(k+1) + edges(k))/2;
        end

        subplot(2,nchains,c+nchains);
        plot( bin_center, avg_energy, '-b', 'linewidth', 2 );
        axis([param_range(1,p), param_range(2,p), -Inf, Inf ]);
        xlabel( sprintf('%s', cfg.param_names{p}) );
        ylabel( 'mean energy' );

    end

    for c = 1:nchains
        subplot(2,nchains,c);
        axis( [param_range(1,p), param_range(2,p), 0, maxpdf]);
    end

    pause;
    
end
