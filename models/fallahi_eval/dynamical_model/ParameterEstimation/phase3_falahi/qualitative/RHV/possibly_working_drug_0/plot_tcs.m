% tc(:,1) = log2(possiblyworkingdrug000001(:,2)./possiblyworkingdrug000001(1,2));
% tc(:,2) = log2(possiblyworkingdrug000002(:,2)./possiblyworkingdrug000002(1,2));
% tc(:,3) = log2(possiblyworkingdrug000003(:,2)./possiblyworkingdrug000003(1,2));
% tc(:,4) = log2(possiblyworkingdrug000004(:,2)./possiblyworkingdrug000004(1,2));
% tc(:,5) = log2(possiblyworkingdrug000005(:,2)./possiblyworkingdrug000005(1,2));
% tc(:,6) = log2(possiblyworkingdrug000006(:,2)./possiblyworkingdrug000006(1,2));
% tc(:,7) = log2(possiblyworkingdrug000007(:,2)./possiblyworkingdrug000007(1,2));
% tc(:,8) = log2(possiblyworkingdrug000008(:,2)./possiblyworkingdrug000008(1,2));
% tc_(:,9) = log2(possiblyworkingdrug000009(:,2)./possiblyworkingdrug000009(1,2));
% tc(:,10) = log2(possiblyworkingdrug000010(:,2)./possiblyworkingdrug000010(1,2));
% tc(:,11) = log2(possiblyworkingdrug000011(:,2)./possiblyworkingdrug000011(1,2));
% tc(:,12) = log2(possiblyworkingdrug000012(:,2)./possiblyworkingdrug000012(1,2));
% tc(:,13) = log2(possiblyworkingdrug000013(:,2)./possiblyworkingdrug000013(1,2));
% tc(:,14) = log2(possiblyworkingdrug000014(:,2)./possiblyworkingdrug000014(1,2));


load('lox_results.mat')
figure
 hold on;
 for i = 1:14
    plot(1:length(tc(1:end,1)),tc(1:end,i),'LineWidth',2)
    hold on
end
        % Specify graphics
        xlim([-.01 40]);
        ylim([-1 3 ]);
%         title(tit{j},'fontsize',28);
%         ylabel(tit{j},'fontsize',24,'fontweight','b','fontname','Arial');
        ylabel('Log2 Fold Change','fontsize',24,'fontweight','b','fontname','Arial');
        xlabel('Time','fontsize',24,'fontweight','b','fontname','Arial');
        set(gca,'fontsize',20,'fontname','Arial');
        set(gca,'xticklabel',[])

        set(gcf,'color','w');
        box on;
%         set(gca,'XTick',0:20:Tf);



clear
load('rvh_results.mat')
figure
 hold on;
 for i = 1:14
    plot(1:length(tc(1:end,1)),tc(1:end,i),'LineWidth',2)
    hold on
end
        % Specify graphics
        xlim([-.01 400]);
        ylim([-5 3 ]);
%         title(tit{j},'fontsize',28);
%         ylabel(tit{j},'fontsize',24,'fontweight','b','fontname','Arial');
        ylabel('Log2 Fold Change','fontsize',24,'fontweight','b','fontname','Arial');
        xlabel('Time','fontsize',24,'fontweight','b','fontname','Arial');
        set(gca,'xticklabel',[])

        set(gca,'fontsize',20,'fontname','Arial');
        set(gcf,'color','w');
        box on;
%         set(gca,'XTick',0:20:Tf);