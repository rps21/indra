function [training] = quartile_plots_bigmech(pars,cfg,tag,suppress_plot)

parameters = 10.^(pars');

warning off


skip = 1; 

numtraj = size(parameters,1);

k=size(parameters,1);
j=floor(k/skip);

% each q represents 1 percentile for the area plots 
q1=round(.05*j);
q2=round(.25*j);
q3=round(.5*j);
q4=round(.75*j);
q5=round(.95*j);

% adjust these for the # equations (eqn), end time of simulation (Tf) and #
% of time points you want to simulate (tend)
tend=11; Tf = 172800;  
eqn=8;
tspan = linspace(0,Tf,tend);

% the log_set are the variables we present on log scale
log_set=[];
% dat_set1 is presented log scale, while dat_set2 is presented linearly
dat_set1=[]; dat_set2 = [1,2,3,4,5];

% trajectories = zeros(numtraj,tend,eqn);
    t1 = [cfg.sim_tstart : cfg.sim_dt : cfg.sim_tstop]';   

trajectories = zeros(numtraj,length(t1),eqn); 

n_err = 1;
errindex=0;
for i = 1:numtraj
    par = parameters(i,:);
    i
    % adjust all this solving stuff for your model

    t = [cfg.sim_tstart : cfg.sim_dt : cfg.sim_tstop]';
    tequil = linspace(0,1e5,250)';


    
%       t = linspace(0,1e5,250)';
    init = zeros(1,725);
    init(1) = par(72);
    init(2) = par(55);    %rptor
    init(3) = par(56);    %raf1 count
    init(4) = par(57);    %grb2 count
    init(5) = par(58);    %mapk1 count
    init(6) = par(59);    %sos1 count
    init(7) = par(60);    %mtor count
    init(8) = par(61);    %rps6kb1 count
    init(9) = par(62);    %mapk3 mRNA count
    init(10) = par(63);    %braf  count
    init(11) = par(64);    %akt count
    init(12) = par(65);    %egfr mRNA count
    init(13) = par(66);    %kras count
    init(14) = par(67);    %rictor count
    init(15) = par(68);    %map2k1 count
    init(16) = par(69);    %cdkn1b count    
        init(17) = par(70);    %s6 count    
    init(18) = 0;    %drug count    

       [err, species, obsv_equil] = model_bigmech_test(tequil,init,par);


    clear init

    init = species(end,:);
    if tag == 1
        init(18) = par(71);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         init(2) = par(2);
%     elseif tag == 2
%         init(1) = par(1)*4;
%         init(2) = par(2);
    end

    
    [err, species, obsv] = model_bigmech_test(t,init,par);
    
    index = find(obsv<1e-5);
    obsv(index) = 1e-5; 
    
    obsv(:,1) = log2(obsv(:,1)./obsv(1,1));         
    obsv(:,2) = log2(obsv(:,2)./obsv(1,2));         
    obsv(:,3) = log2(obsv(:,3)./obsv(1,3));         
    obsv(:,4) = log2(obsv(:,4)./obsv(1,4));         
    obsv(:,5) = log2(obsv(:,5)./obsv(1,5));         
    obsv(:,6) = log2(obsv(:,6)./obsv(1,6));         
    obsv(:,7) = log2(obsv(:,7)./obsv(1,7));         
    obsv(:,8) = log2(obsv(:,8)./obsv(1,8));     
    
    if err == 0 
        obsv=obsv';
    	obsv_equil = obsv_equil';
        for m=1:eqn
            bin=contain(m,log_set);
            if (bin)
                  trajectories(i,:,m)=log10(obsv(m,1:end));
            else
%                 size(trajectories)
%                 size(trajectories(i,:,m))
%                 size(obsv)
%                 size(obsv(m,1:end))
                trajectories(i,:,m)=obsv(m,1:end);
                trajectories_equil(i,:,m)=obsv_equil(m,1:end);
            end
        end
    else
        errindex(n_err) = i;
        n_err = n_err+1;
    end
    
end

% sorttraj = sort(trajectories);
trajectories = sort(trajectories);
trajectories_equil = sort(trajectories_equil);

 

% t1=tspan;
t1 = [cfg.sim_tstart : cfg.sim_dt : cfg.sim_tstop]';%t;
tit={'pmek','perk','Akt_p308','Akt_p473','pmtor','ps6k','ps6','pp27'};
% 1  struct('name','p-mek', 'units','log10'), ...
% 2  struct('name','p-erk', 'units','log10'), ...
% 3  struct('name','akt_p308', 'units','log10'), ...
% 4  struct('name','akt_p473', 'units','log10'), ...
% 5  struct('name','p-mtor', 'units','log10'), ...
% 6  struct('name','p-s6k', 'units','log10'), ...
% 7  struct('name','p-s6', 'units','log10'), ...
% 8  struct('name','p-p27', 'units','log10'), ...

% insert a line loading your raw data for plotting, if desired
% Low Dose
if tag == 1
    time_data = cfg.timepoints;
    data  = cfg.data{1}.mean; 
%     data(2,1) = NaN;
    stdev = cfg.data{1}.stdev;
end

% High Dose
% if tag == 2
%     time_data = cfg.timepoints;
%     data  = cfg.data{2}.mean;
%     stdev = cfg.data{2}.stdev;
% end

% bounds on y-axis of each subplot
% lb = zeros(1,10);
lb = -2*ones(1,10);

% lb=[0 0 0 0 0 0 0 0];   

ub=2.1*ones(1,10);
% ub(1) = 2.1;
% ub = [2.1 1.1 1.1 1.1 1.1 1.1 1.1 1.1];

% outercolor = [0.95 0.95 0.95]; % light gray
% innercolor = [0.8 0.8 0.8]; % dark gray

outercolor = [0.9 0.9 1]; % light blue
innercolor = [0 0 0.85]; % dark blue

% if suppress_plot ~= 1
% %     figure;    
% 	j = 1;
%      for m = 1:eqn;
% %         subplot(1,6,j)
%         figure;
%         hold on;
%         % Specify graphics
%         xlim([-.5 (Tf)]);
%         ylim([lb(m) ub(m)]);
% %         title(tit{j},'fontsize',28);
%         ylabel(tit{j},'fontsize',24,'fontweight','b','fontname','Arial');
%         xlabel('Time (seconds)','fontsize',24,'fontweight','b','fontname','Arial');
%         set(gca,'fontsize',20,'fontname','Arial');
%         set(gcf,'color','w');
%         box on;
%         
% %         set(gca,'XTick',0:20:Tf);
%         % plot ensemble
% %         plot(t1,trajectories(1800:20:q5,:,m),'-r','LineWidth',1); 
%         area(tequil,trajectories_equil(q5,:,m),'FaceColor',outercolor,'LineStyle','none'); hold on;        
%         area(tequil,trajectories_equil(q4,:,m),'FaceColor',innercolor,'LineStyle','none'); hold on;       
%         area(tequil,trajectories_equil(q2,:,m),'FaceColor',outercolor,'LineStyle','none'); hold on;
%         area(tequil,trajectories_equil(q1,:,m),'FaceColor','w','LineStyle','none'); hold on;
%         plot(tequil,trajectories_equil(q3,:,m),'-r','LineWidth',3);       
%         % plot data - comment out if not plotting raw data
%         bin=contain(m,dat_set1); % log scale data
%         if (bin)
%             index=contfind(m,dat_set1); 
%             % change variable names here to reflect the data you loaded
%             errorbar(time_data,data(:,1),stdev(:,1),'+k','LineWidth',3.0);
%         end
%         bin=contain(m,dat_set2); % linear scale data
%         if (bin)
%             index=contfind(m,dat_set2);    
%             errorbar(time_data,data(:,m),stdev(:,m),'+k','LineWidth',3.0);
% %            end
%         end
%         hold off;
% 	j = j+1;
%     end
% 
% end


if suppress_plot ~= 1
%     figure;    
	j = 1;
     for m = 1:7;
%         subplot(1,6,j)
        figure;
        hold on;
        % Specify graphics
        xlim([-.01 (Tf)./(60*60)]);
        ylim([lb(m) ub(m)]);
%         title(tit{j},'fontsize',28);
        ylabel(tit{j},'fontsize',24,'fontweight','b','fontname','Arial');
        xlabel('Time (hours)','fontsize',24,'fontweight','b','fontname','Arial');
        set(gca,'fontsize',20,'fontname','Arial');
        set(gcf,'color','w');
        box on;
%         set(gca,'XTick',0:20:Tf);
        % plot ensemble
%         plot(t1,trajectories(1800:20:q5,:,m),'-r','LineWidth',1); 

        area(t1./(60*60),trajectories(q5,:,m),'FaceColor',outercolor,'LineStyle','none'); hold on;        
        area(t1./(60*60),trajectories(q4,:,m),'FaceColor',innercolor,'LineStyle','none'); hold on;       
        area(t1./(60*60),trajectories(q2,:,m),'FaceColor',outercolor,'LineStyle','none'); hold on;
        area(t1./(60*60),trajectories(q1,:,m),'FaceColor','w','LineStyle','none'); hold on;
        plot(t1./(60*60),trajectories(q3,:,m),'-r','LineWidth',3);       
        % plot data - comment out if not plotting raw data
        bin=contain(m,dat_set1); % log scale data
        if (bin)
            index=contfind(m,dat_set1); 
            % change variable names here to reflect the data you loaded
            errorbar(time_data./(60*60),data(:,1),stdev(:,1),'+k','LineWidth',3.0);
        end
        bin=contain(m,dat_set2); % linear scale data
        if (bin)
            index=contfind(m,dat_set2);    
            errorbar(time_data./(60*60),data(:,m),stdev(:,m),'+k','marker','s','markersize',9,'markerfacecolor','black','LineWidth',3.0);
%            end
        end
        hold off;
	j = j+1;
    end

end
end

function bin=contain(m,set)
    bin=0;
    k=size(set,2);
    for j=1:k
        if (m==set(j))
            bin=1;    
            return;
        end
    end
    return;
end
    
function index=contfind(m,set)
index=0;
k=size(set,2);    
    for j=1:k
        if (m==set(j))
            index=j;    
            return;
        end
    end
    return;
end
