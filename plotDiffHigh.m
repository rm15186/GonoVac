 figure('name','Comparison of Strategies');
 hold on;
 days = [-500:4000];
    plot(days,plot_change_either1High,'color',[255/255,128/255,0],'linewidth',1);
    %plot(days,plot_change1(:,1),'b-');
    %plot(days,plot_change1(:,2),'r-');
    plot(days,plot_change_either2High,'linewidth',1,'color',[19/255,145/255,62/255]);
    %plot(days,plot_change2(:,1),'b--');
    %plot(days,plot_change2(:,2),'r--');
    plot(days,plot_change_either3High,'color',[0.494,0.184,0.556],'linewidth',1);
    xline(0);
    %plot(days,plot_change3(:,1),'b:');
    %plot(days,plot_change3(:,2),'r:');
    xlabel('Time (days)','fontsize',14,'Interpreter','latex');
    ylabel('Prevalence (%)','fontsize',14,'Interpreter','latex');
    legend('Strategy 1 - Vaccination in Childhood','Strategy 2 - Vaccination at Screening','Strategy 3 - Vaccination at Diagnosis','fontsize',14,'Interpreter','latex');
    
    box on;
    grid on;
    
    figure('name','comparison of vaccinated percentages');
    hold on;
        %avg_vac_current3 = 100*avg_vac_current/N;
        %std_vac_current3 = 100^2*std_vac_current/N^2; %normal
        blue = [0,0.447,0.7410];
        green = [19/255,145/255,62/255];
        purple = [0.494,0.184,0.556];
        plot([0:n_Days],avg_vac_current1High,'color',[255/255,128/255,0],'linewidth',1);
        plot([0:n_Days],avg_vac_current2High,'linewidth',1,'color',[19/255,145/255,62/255]);
        plot([0:n_Days],avg_vac_current3High,'color',[0.494,0.184,0.556],'linewidth',1);
        %shadedErrorBar([0:n_Days],avg_vac_current1,[std_vac_current1(:,1),std_vac_current1(:,1)]);
        %shadedErrorBar([0:n_Days],avg_vac_current2,[std_vac_current2(:,1),std_vac_current2(:,1)]);
        %shadedErrorBar([0:n_Days],avg_vac_current3,[std_vac_current3(:,1),std_vac_current3(:,1)]);
        xlabel('Time (days)','fontsize',14,'Interpreter','latex');
        ylabel('Average percentage of people protected','fontsize',14,'Interpreter','latex');
        %title('average no of people protected','fontsize',14,'Interpreter','latex');
        legend('Strategy 1 - Vaccination in Childhood','Strategy 2 - Vaccination at Screening','Strategy 3 - Vaccination at Diagnosis','fontsize',14,'Interpreter','latex');
        box on;
    grid on;
    
    figure('name','comparison of drug doses');
    hold on;
        plot((0:n_Days),avg_cefta1High,'color',[255/255,128/255,0],'linewidth',1);
        plot((0:n_Days),avg_cefta2High,'color',[19/255,145/255,62/255],'linewidth',1);
        plot((0:n_Days),avg_cefta3High,'color',[0.494,0.184,0.556],'linewidth',1);
        xlabel('Time (days)','fontsize',14,'Interpreter','latex');
        ylabel('Cumulative drug doses administered','fontsize',14,'Interpreter','latex');
        legend('Strategy 1 - Vaccination in Childhood','Strategy 2 - Vaccination at Screening','Strategy 3 - Vaccination at Diagnosis','fontsize',14,'Interpreter','latex');
        box on;
        grid on;
        
   figure('name','comparison of vaccination doses');
   hold on;
        plot((0:n_Days),avg_vac_doses1High,'color',[255/255,128/255,0],'linewidth',1);
        plot((0:n_Days),avg_vac_doses2High,'color',[19/255,145/255,62/255],'linewidth',1);
        plot((0:n_Days),avg_vac_doses3High,'color',[0.494,0.184,0.556],'linewidth',1);
        xlabel('Time (days)','fontsize',14,'Interpreter','latex');
        ylabel('Cumulative vaccine doses administered','fontsize',14,'Interpreter','latex');
        legend('Strategy 1 - Vaccination in Childhood','Strategy 2 - Vaccination at Screening','Strategy 3 - Vaccination at Diagnosis','fontsize',14,'Interpreter','latex');
        box on;
        grid on;
        
        
        
        