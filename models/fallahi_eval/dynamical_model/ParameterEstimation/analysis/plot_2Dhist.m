function [] = plot_2Dhist( logsamples, cfg );
%PLOT_2DHIST plot 2-d histograms of samples
%
%  [] = plot_2Dhist( logsamples, cfg )
%  
%  where 'logsamples' is a (N x P) array of samples,
%  with N=number of samples, P=number of parameters; and
%  'cfg' is the configuration struct.
%
%  Histograms are shown for one parameter at a time. Hit any key to see
%  the next histogram.

param_range = cfg.param_range;

nsamples = size(logsamples,1);
nparams = size(logsamples,2);

res = 16;
vplots = 4;
hplots = ceil(nparams/vplots);

for p1 = 1:nparams

    % skip params with singleton range
    if param_range(1,p1) == param_range(2,p1)
        continue;
    end

    fh = figure(1);
    set( fh, 'Color', [1 1 1] );
    %title( sprintf('%s', cfg.param_names{p1} ), 'fontsize', 11 );

    for p2 = 1:nparams

        sph = subplot(vplots,hplots,p2);
        cla(sph);

        % skip self
        if p1 == p2
            continue;
        end
 
        % skip params with singleton range
        if param_range(1,p2) == param_range(2,p2)
            continue;
        end

        
        binedges{1} = linspace( param_range(1,p2), param_range(2,p2), res+1 );
        binedges{2} = linspace( param_range(1,p1), param_range(2,p1), res+1 );

        [N, C] = hist3( [logsamples(:,p2), logsamples(:,p1)], 'Edges', binedges );
        C{1} = C{1}(1:end-1);
        C{2} = C{2}(1:end-1);
        N = N(1:end-1,1:end-1);

        pdf = N/( (binedges{1}(2)-binedges{1}(1))*(binedges{2}(2)-binedges{2}(1))*sum(sum(N)) );
        maxpdf = max(max(pdf));

        colorscale = linspace(0,maxpdf,20);
        [~,ch] = contourf( C{1}, C{2}, pdf', colorscale );
        caxis([0,maxpdf]);
        set(ch,'edgecolor','none');
        xlabel( sprintf('%s', cfg.param_names{p2} ), 'fontsize', 9 );

        if (mod(p2,hplots)==1)
            ylabel( sprintf('%s', cfg.param_names{p1} ), 'fontsize', 9 );
        end

    end
    
    pause;
    
end
