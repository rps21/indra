function [] = plot_swap( swap_acceptance, beta_history, cfg, fig_idx )
%PLOT_SWAP
%
%  [] = plot_swap( swap_acceptance, beta_history, cfg, fig_idx )

n_chains = size( swap_acceptance, 1 )+1;
n_adjust = floor( size(swap_acceptance,2) / cfg.adapt_beta_interval );

steps = [1:n_adjust];

avg_swap_acceptance = zeros( n_chains-1, n_adjust );
for i = 1: n_adjust
    avg_swap_acceptance(:,i) = mean( swap_acceptance(:,(i-1)*cfg.adapt_beta_interval+1:i*cfg.adapt_beta_interval), 2);
end
avg_swap_acceptance = avg_swap_acceptance';

temperature = zeros(n_chains, n_adjust);
for i = 1: n_adjust
    temperature(:,i) = 1./(beta_history(:, (i-1)*cfg.adapt_beta_interval+1));
end
temperature = temperature';


fh = figure(fig_idx);
set( fh, 'Color', [1 1 1] );
for i = 1 : (n_chains-1)
    subplot( 2, n_chains-1, i );
    stairs( steps, temperature(:,i+1),  '-r' );
    ylabel( 'temperature', 'fontSize',11 );
    set(gca,'yscale','log');
    axis( [0, n_adjust, 1, 20] );

    subplot( 2, n_chains-1, i+(n_chains-1) );
    stairs( steps, avg_swap_acceptance(:,i), '-b' );
    ylabel( 'swap accept', 'fontSize',11 );
    set(gca,'yscale','log');
    axis( [0 n_adjust 0.001, 1] );
    hold on;
    plot( steps, 0.22*ones(size(steps)), '--k' );
    hold off;

%    subplot( 4, n_chains-1, i+(n_chains-1)*2 );
%    stairs( steps(1:end-1), log2(avg_swap_acceptance(2:end,i)./avg_swap_acceptance(1:end-1,i)),  '-b' );
%    ylabel( 'log(swap accept ratio)', 'fontSize',11 );
%    hold on;
%    plot( steps, zeros(size(steps)), '--k' );
%    axis( [0 n_adjust -4, 4] );
%    hold off;

%    subplot( 4, n_chains-1, i+(n_chains-1)*3 );
%    stairs( steps(1:end-1), log2(temperature(2:end,i+1)./temperature(1:end-1,i+1)),  '-r' );
%    xlabel( sprintf( 'step\n chain %d, swap %d,%d', i+1, i, i+1 ), 'fontSize', 12 );
%    ylabel( 'log(temp ratio)', 'fontSize',11 );
%    hold on;
%    plot( steps, zeros(size(steps)), '--k' );
%    axis( [0 n_adjust -0.5, 0.5] );
%    hold off;

end


fh = figure(fig_idx+1);
set( fh, 'Color', [1 1 1] );
colors = {'r','y','g','b','c'};
linestyle = {'-r','-m','-g','-b','-c'};
for i = 1:(n_chains-1)

    logT   = log2(temperature(2:end,i+1)./temperature(1:end-1,i+1)) ...
             - log2(temperature(2:end,i)./temperature(1:end-1,i));
    
    logSwA = -cfg.adapt_beta_rate*log2(avg_swap_acceptance(2:end,i)./avg_swap_acceptance(1:end-1,i));

    goodidx = find( and( ~isnan(logSwA), ~isinf(logSwA) ) );
    logT   = logT(goodidx);
    logSwA = logSwA(goodidx);

    scatter( logT, logSwA, 9, colors{i} );
    hold on;

    coeff = polyfit( logT, logSwA, 1 );
    x = [-2:.01:2];
    y = coeff(1)*x + coeff(2);
    plot(x, y, linestyle{i});

end
plot( x,zeros(size(x)), '--k' );
hold off;
axis( [-0.5 0.5, -0.6 0.6] );
xlabel( 'log( T_{n+1} / T_n )', 'fontsize',11 );
ylabel( '-\alpha log( SwA_{n+1} / SwA_n )', 'fontsize',11 );





