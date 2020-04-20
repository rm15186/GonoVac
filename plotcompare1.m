figure('name','Comparison of Strategies');
 hold on;
 days = [-500:4000];
    plot(days,plot_change_either1,'color',[255/255,128/255,0],'linewidth',1);
    %plot(days,plot_change1(:,1),'b-');
    %plot(days,plot_change1(:,2),'r-');
    plot(days,plot_change_either2,'linewidth',1,'color',[19/255,145/255,62/255]);
    %plot(days,plot_change2(:,1),'b--');
    %plot(days,plot_change2(:,2),'r--');
    plot(days,plot_change_either3,'color',[0.494,0.184,0.556],'linewidth',1);
    xline(0);
    %plot(days,plot_change3(:,1),'b:');
    %plot(days,plot_change3(:,2),'r:');
    xlabel('Time (days)','fontsize',14,'Interpreter','latex');
    ylabel('Prevalence (%)','fontsize',14,'Interpreter','latex');
    legend('Strategy 1 - Vaccination in Childhood','Strategy 2 - Vaccination at Screening','Strategy 3 - Vaccination at Diagnosis','fontsize',14,'Interpreter','latex');
    
    box on;
    grid on;