function [] = plot_diffhist( samples1, samples2, cfg  );
%plot_DIFFHIST plot difference in 2-d histograms for parameter samples
%
%  [] = plot_DIFFHIST( samples1, samples2, cfg )
%  
%  where 'samples1' is a (S x P) array of samples,
%  'samples2' is another (S x P) array of samples,
%  with S=number of samples, P=number of parameters; and
%  'cfg' is the configuration struct.
%
%  Histograms are shown for one parameter at a time. Hit any key to see
%  the next histogram.
width = 1.6;
param_range = [ cfg.param_location - width*cfg.param_scale; cfg.param_location + width*cfg.param_scale ];
%param_range = cfg.param_range;

nsamples = size(samples1,1);
nparams = size(samples1,2);

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

        [N1, C] = hist3( [samples1(:,p2), samples1(:,p1)], 'Edges', binedges );
        [N2, ~] = hist3( [samples2(:,p2), samples2(:,p1)], 'Edges', binedges );

        C{1} = C{1}(1:end-1);
        C{2} = C{2}(1:end-1);

        N1 = N1(1:end-1,1:end-1);
        pdf1 = N1 / sum(sum(N1));

        N2 = N2(1:end-1,1:end-1);
        pdf2 = N2 / sum(sum(N2));
        
        maxpdf1 = max(max(pdf1));
        maxpdf2 = max(max(pdf2));

        maxpdf = max( [maxpdf1, maxpdf2] );

        colorscale = linspace(-1,1,20);
        [~,ch] = contourf( C{1}, C{2}, ((pdf2-pdf1)/maxpdf)', colorscale );
        caxis([-1,1]);
        set(ch,'edgecolor','none');
        xlabel( sprintf('%s', cfg.param_names{p2} ), 'fontsize', 9 );

        if (mod(p2,hplots)==1)
            ylabel( sprintf('%s', cfg.param_names{p1} ), 'fontsize', 9 );
        end

    end
    
    pause;
    
end
