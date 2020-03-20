% plot_processed data in the normal style with std for everything except
% prevalence post vaccine
N = 10000;
n_Days = 4000;

            close all
            %plot 25th and 75th percentiles, not standard deviation as we
            %arent going to have negative prevalence so normal distribution
            %probably doesnt make sense here
            figure('name', 'Average Prevalence with interquartile range');
                hold on;
                
                shadedErrorBar([0:n_Days],plot_either,[plot_either-i25e,i75e-plot_either],'lineprops','k');
                shadedErrorBar([0:n_Days],plot_data(:,1),[plot_data(:,1)-i25(:,1),i75(:,1)-plot_data(:,1)],'lineprops','b');
                shadedErrorBar([0:n_Days],plot_data(:,2),[plot_data(:,2)-i25(:,2),i75(:,2)-plot_data(:,2)],'lineprops','r');
                legend('Either Strain','Non-AMR','AMR','fontsize',14,'Interpreter','latex');
                xlabel('Time (days)','fontsize',14,'Interpreter','latex')
                ylabel('Average Prevalence (%)','fontsize',14,'Interpreter','latex');
                box on;
                grid on;
            
            %add the overall prevalence and make these thicker
            figure('name', 'Average Prevalence');
                hold on;
                plot([0:n_Days],plot_data(:,1),'b-');
                plot([0:n_Days],plot_data(:,2),'r-');
                plot([0:n_Days],plot_either,'k-');
                legend('non-AMR','AMR','Total Prevalence','fontsize',14,'Interpreter','latex');
                xlabel('Time (days)','fontsize',14,'Interpreter','latex')
                ylabel('Average Prevalence (%)','fontsize',14,'Interpreter','latex');
                box on;
                grid on;
            
            
            
            % drug administration of each drug given by the cumulative sum
            % of the daily dosage of each drug
                figure('name','Dosage','color','w');
                    hold on;
                   
                    %shadedErrorBar([0:n_Days], cumsum(avg_cipr_doses(:,1)),[std_cipr_doses(:,:,1),std_cipr_doses(:,:,1)]);
                    shadedErrorBar([0:n_Days],cumsum(avg_cefta_doses(:,1)),[std_cefta_doses(:,:,1),std_cefta_doses(:,:,1)]);
                    %legend('Cipr/A','Ceft/A','location','northwest','fontsize',14,'Interpreter','latex');
                    xlabel('Time (days)','fontsize',14,'Interpreter','latex');
                    ylabel('Number of doses','fontsize',14,'Interpreter','latex')
                    title('Cumulative drug doses administered','fontsize',14,'Interpreter','latex');
                    box on;
                    grid on;
                    
  
                
                figure('name','Average doses of vaccine');
                    hold on;
                    %plot([0:n_Days], cumsum(avg_vac_doses));
                    shadedErrorBar([0:n_Days],cumsum(avg_vac_doses),[std_vac_doses(:,1),std_vac_doses(:,1)])
                    xlabel('Time (days)','fontsize',14,'Interpreter','latex');
                    ylabel('No of vaccine doses given','fontsize',14,'Interpreter','latex');
                    title('average vaccine doeses given','fontsize',14,'Interpreter','latex');
                    box on;
                    grid on;
                    
                figure('name','average no. of vaccinated people');
                    hold on;
                    avg_vac_current = 100*avg_vac_current/N;
                    std_vac_current = 100^2*std_vac_current/N^2; %normalise
                    plot([0:n_Days], avg_vac_current);
                    shadedErrorBar([0:n_Days],avg_vac_current,[std_vac_current(:,1),std_vac_current(:,1)])
                    xlabel('Time (days)','fontsize',14,'Interpreter','latex');
                    ylabel('Average percentage of people protected','fontsize',14,'Interpreter','latex');
                    %title('average no of people protected','fontsize',14,'Interpreter','latex');
                    box on;
                grid on;
                
                %TODO plot end of burn in and then plot effect of vaccine
                %currently not working but should just plot what happens
                %over the 2000 day burn in 
                 figure('name','burn in');
                     hold on;
                     %plot([0:3000],plot_burn_in_prev_either,'k');
                     %plot([0:3000],plot_burn_in_prev(:,1),'b');
                     %plot([0:3000],plot_burn_in_prev(:,2),'r');
                     size([0:3000]);
                     size(std_burn_in_prev(:,1));
                     size(plot_burn_in_prev);
                     size(std_burn_in_prev_either);
                     shadedErrorBar([0:3000],plot_burn_in_prev_either,[std_burn_in_prev_either],'lineprops','k')
                     shadedErrorBar([0:3000],plot_burn_in_prev(:,1),[std_burn_in_prev(:,1),std_burn_in_prev(:,1)],'lineprops','b')
                     shadedErrorBar([0:3000],plot_burn_in_prev(:,2),[std_burn_in_prev(:,2),std_burn_in_prev(:,2)],'lineprops','r')
                     xlabel('Time (days)','fontsize',14,'Interpreter','latex');
                     ylabel('prevalence in burn in','fontsize',14,'Interpreter','latex');
                     title('prevalence over burn in period','fontsize',14,'Interpreter','latex');
                     box on;
                     grid on;
                     
                 
                 %stitch the data together 
                 days = [0:500+n_Days]; %end of burn in then length of simulation
                 burn = plot_burn_in_prev_either([2500:2999]);
                 burn_strains = plot_burn_in_prev([2500:2999],:);
                 %size(burn_strains)
                 %size(plot_data)
                 plot_change = cat(1,burn_strains,plot_data);
                 plot_change_either = cat(1,burn,plot_either);
                 size(plot_change);
                 size(days);
                 figure('name','steady state then vaccine');
                    hold on;
                    plot(days,plot_change_either,'k-');
                    plot(days,plot_change(:,1),'b-');
                    plot(days,plot_change(:,2),'r-')
                    xline(500);
                    xlabel('Time (days)','fontsize',14,'Interpreter','latex');
                    ylabel('Prevalence (%)','fontsize',14,'Interpreter','latex');
                    title('Impact of vaccine on prevalence','fontsize',14,'Interpreter','latex');
                    legend('Eiter Strain','Non-AMR','AMR','fontsize',14,'Interpreter','latex')
                    box on;
                    grid on;
                 