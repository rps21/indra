function [] = plot_relstep( step_acceptance, relstep_history, cfg, fig_idx )
%PLOT_RELSTEP
%
%  [] = plot_relstep( step_acceptance, relstep_history, cfg, fig_idx )


n_chains = size( step_acceptance, 1 );
n_adjust = floor( size(step_acceptance,2) / cfg.adapt_relstep_interval );



steps = [1:n_adjust];

avg_step_acceptance = zeros( n_chains, n_adjust );
for i = 1: n_adjust
    avg_step_acceptance(:,i) = mean( step_acceptance(:,(i-1)*cfg.adapt_relstep_interval+1:i*cfg.adapt_relstep_interval), 2) / cfg.nsteps;
end
avg_step_acceptance = avg_step_acceptance';

relstep = zeros(n_chains, n_adjust);
for i = 1: n_adjust
    relstep(:,i) = relstep_history(:, (i-1)*cfg.adapt_relstep_interval+1);
end
relstep = relstep';


fh = figure(fig_idx);
set( fh, 'Color', [1 1 1] );
for i = 1 : (n_chains)
    subplot( 2, n_chains, i );
    stairs( steps, relstep(:,i),  '-r' );
    ylabel( 'rel step', 'fontSize',11 );
    set(gca,'yscale','log');
    axis( [0, n_adjust, 1e-4, 1e-1] );

    subplot( 2, n_chains, i+n_chains );
    stairs( steps, avg_step_acceptance(:,i), '-b' );
    ylabel( 'step accept', 'fontSize',11 );
    set(gca,'yscale','log');
    axis( [0, n_adjust, 0.005, 1] );
    hold on;
    plot( steps, 0.22*ones(size(steps)), '--k' );
    hold off;

%    subplot( 4, n_chains, i+n_chains*2 );
%    stairs( steps(1:end-1), log2(avg_step_acceptance(2:end,i)./avg_step_acceptance(1:end-1,i)),  '-b' );
%    ylabel( 'log(step accept ratio)', 'fontSize',11 );
%    hold on;
%    plot( steps, zeros(size(steps)), '--k' );
%    axis( [0 n_adjust -4, 3] );
%    hold off;

%    subplot( 4, n_chains, i+n_chains*3 );
%    stairs( steps(1:end-1), log2(relstep(2:end,i)./relstep(1:end-1,i)),  '-r' );
%    xlabel( sprintf( 'step\n chain %d', i ), 'fontSize', 12 );
%    ylabel( 'log(rel step ratio)', 'fontSize',11 );
%    hold on;
%    plot( steps, zeros(size(steps)), '--k' );
%    axis( [0 n_adjust -1, 1] );
%    hold off;

end


fh = figure(fig_idx+1);
set( fh, 'Color', [1 1 1] );
colors = {'r','y','g','b','c'};
linestyle = {'-r','-m','-g','-b','-c'};
for i = 1 : n_chains

    logRS  = log2(relstep(2:end,i)./relstep(1:end-1,i));
    logStA = -cfg.adapt_relstep_rate*log2(avg_step_acceptance(2:end,i)./avg_step_acceptance(1:end-1,i));

    scatter( logRS, logStA, 9, colors{i} );
    hold on;

    coeff = polyfit( logRS, logStA, 1 );
    x = [-2:.01:2];
    y = coeff(1)*x + coeff(2);
    plot(x, y, linestyle{i});

end
plot( x, zeros(size(x)), '--k' );
hold off;
axis( [-1, 1, -1, 1] );
xlabel( 'log( rs_{n+1} / rs_n )' );
ylabel( '-\alpha log( StA_{n+1} / StA_n )' );



