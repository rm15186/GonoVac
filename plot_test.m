



% plot_data = mean(all_data,3);
%             plot_either = mean(all_either,3);
%             
% %             size(all_burn_in_prev)
% %             size(all_burn_in_prev_either)
%             plot_burn_in_prev_either = mean(all_burn_in_prev_either,3); %2001?
%             plot_burn_in_prev = mean(all_burn_in_prev,3);
%             
%             range_burn_prev = [quantile(all_burn_in_prev,0.25,2),quantile(all_burn_in_prev,0.75,2)];
%             range_burn_prev_either = [quantile(all_burn_in_prev_either,0.25,2),quantile(all_burn_in_prev_either,0.75,2)];
%             
%             %size(range_burn_prev)
%             %size(range_burn_prev_either)
% %             size(range_burn_either)
%             
%             
%              range_burn_either = zeros(3001,2);
% %             range_burn_in_prev = [quantile(all_burn_in_prev,0.25,2),quantile(all_burn_in_prev,0.75,2)];
% %             range_burn_in_prev_either = [quantile(all_burn_in_prev_either,0.25,3),quantile(all_burn_in_prev_either,0.75,3)];
%              plot_range_burn_either(:,1) = plot_burn_in_prev_either-range_burn_prev_either(:,1);
%              plot_range_burn_either(:,2) = range_burn_in_prev_either(:,2)-plot_burn_in_prev_either;
%              plot_range_burn_either(:,1) = range_burn_prev_either(:,1);
%              plot_range_burn_either(:,2) = range_burn_in_prev_either(:,2);
%              
%              plot_range_burn_strains = [plot_burn_in_prev_either-range_burn_prev(:,1),range_burn_prev(:,2)-plot_burn_in_prev_either];
% %             
% %              size(plot_burn_in_prev_either)
% %              size(range_burn_in_prev_either(:,1))
% %              size(plot_burn_in_prev_either-range_burn_in_prev_either(:,1))
% %              size(plot_range_burn_either)
% %              size(plot_burn_in_prev)
% %              size(plot_burn_in_prev(:,2))
%             
%             %standard deviation at all points for plotting confidence
%             %intervals, if we want to do that
%             %conf = std(all_data,0,3);
%             %confe = std(all_either,0,2);
%             %size(conf);
%             %quartiles for error bars that dont go below 0
%             
%             %% unsupress these
%             i25 = quantile(all_data,0.25,3); %25th percentile
%             i75 = quantile(all_data,0.75,3);%75th percentile
%             
%             i25e = quantile(all_either,0.25,3); %2 or 3? 3 works!
%             i75e = quantile(all_either,0.75,3);
%   
%             
%             avg_vac_current = mean(all_vac_current,2);
%             %std_vac_current = std(all_vac_current,0,2);
%             range_vac_current = [quantile(all_vac_current,0.25,2),quantile(all_vac_current,0.75,2)];
%             plot_range_vac_current = [avg_vac_current-range_vac_current(:,1),range_vac_current(:,2)-avg_vac_current];
%             
%             avg_vac_doses = mean(all_vac_doses,2);
%             %std_vac_doses = std(all_vac_doses,0,2);
%             range_vac_doses = [quantile(all_vac_doses,0.25,2),quantile(all_vac_doses,0.75,2)];
%             plot_range_vac_doses = [avg_vac_doses-range_vac_doses(:,1),avg_vac_doses-range_vac_doses(:,2)-avg_vac_doses];
%            
%             % avg_cipr_doses = mean(all_cipr_doses,2); %we don't use these
%             %std_cipr_doses = std(all_cipr_doses,0,2);
%             avg_cefta_doses = mean(all_cefta_doses,2);
%             %std_cefta_doses = std(all_cefta_doses,0,2);
%             range_cefta_doses = [quantile(all_cefta_doses,0.25,2),quantile(all_cefta_doses,0.75,2)];
%             plot_range_cefta_doses = [avg_cefta_doses-range_cefta_doses(:,1),range_cefta_doses(:,2)-avg_cefta_doses];

         
            %% plot
            close all 
            figure('name', 'Average Prevalence with interquartile range');
                hold on;
                shadedErrorBar([0:n_Days],plot_either,[plot_either-i25e,i75e-plot_either],'lineprops','k');
                shadedErrorBar([0:n_Days],plot_data(:,1),[plot_data(:,1)-i25(:,1),i75(:,1)-plot_data(:,1),],'lineprops','b');
                shadedErrorBar([0:n_Days],plot_data(:,2),[plot_data(:,2)-i25(:,2),i75(:,2)-plot_data(:,2)],'lineprops','r');
                legend('Either Strain','Non-AMR','AMR','fontsize',14,'Interpreter','latex');
                xlabel('Time (days)','fontsize',14,'Interpreter','latex')
                ylabel('Average Prevalence (%)','fontsize',14,'Interpreter','latex');
                box on;
                grid on;
            
            %add the overall prevalence and make these thicker
            figure('name', 'Average Prevalence');
                hold on;
                plot([0:n_Days],plot_data(:,1),'b-','LineWidth',1);
                plot([0:n_Days],plot_data(:,2),'r-','LineWidth',1);
                plot([0:n_Days],plot_either,'k-','LineWidth',1);
                legend('non-AMR','AMR','Either Strain','fontsize',14,'Interpreter','latex');
                xlabel('Time (days)','fontsize',14,'Interpreter','latex')
                ylabel('Average Prevalence (%)','fontsize',14,'Interpreter','latex');
                box on;
                grid on;
                    
            
            %plot the no. of doses of Cefta, uncomment to plot cipr too
                figure('name','Dosage','color','w');
                    hold on;
                    %shadedErrorBar([0:n_Days], cumsum(avg_cipr_doses(:,1)),[std_cipr_doses(:,:,1),std_cipr_doses(:,:,1)]);
                    shadedErrorBar([0:n_Days],cumsum(avg_cefta_doses(:,1)),plot_range_cefta_doses);
                    %legend('Cipr/A','Ceft/A','location','northwest');
                    xlabel('Time (days)','fontsize',14,'Interpreter','latex');
                    ylabel('Number of doses','fontsize',14,'Interpreter','latex')
                    title('Cumulative Ceftriaxone Doses Administered','fontsize',14,'Interpreter','latex');
                    box on;
                    grid on;
                    
                
                figure('name','Average doses of vaccine');
                    hold on;
                    shadedErrorBar([0:n_Days],cumsum(avg_vac_doses),plot_range_vac_doses);
                    xlabel('Time (days)','fontsize',14,'Interpreter','latex');
                    ylabel('No of vaccine doses given','fontsize',14,'Interpreter','latex');
                    title('Average Vaccine Doses Administered','fontsize',14,'Interpreter','latex');
                    box on;
                    grid on;
                    
                figure('name','Average percentage vaccinated');
                    hold on;
                    shadedErrorBar([0:n_Days],avg_vac_current,plot_range_vac_current);
                    xlabel('Time (days)','fontsize',14,'Interpreter','latex');
                    ylabel('Percentage of people','fontsize',14,'Interpreter','latex');
                    title('Average percentage of people protected by the vaccine','fontsize',14,'Interpreter','latex');
                    box on;
                grid on;
               
                %TODO error bar this isnt quartiles becasue that was messy
                  figure('name','burn in');
                     hold on;
                     %plot([0:3000],plot_burn_in_prev_either,'k');
                     %plot([0:3000],plot_burn_in_prev(:,1),'b');
                     %plot([0:3000],plot_burn_in_prev(:,2),'r');
                     size([0:3000])
                     size(std_burn_in_prev(:,1))
                     size(plot_burn_in_prev)
                     size(std_burn_in_prev_either)
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
                 burn_strains = plot_burn_in_prev([2500:2999],:)';
                 %plot_change = cat(1,burn,plot_data);
                 plot_change = cat(1,burn,plot_either);
                 size(burn_strains);
                 size(plot_data);
                 plot_change_both = cat(1,burn_strains,plot_data);
                 size(plot_change);
                 size(days);
                 figure('name','steady state then vaccine');
                    hold on;
                    plot(days,plot_change,'k');
                    plot(days,plot_change_both(:,1),'b');
                    plot(days,plot_change_both(:,2),'r');
                    xline(500);
                    xlabel('Time (days)','fontsize',14,'Interpreter','latex');
                    ylabel('Prevalence (%)','fontsize',14,'Interpreter','latex');
                    title('Impact of vaccine on prevalence','fontsize',14,'Interpreter','latex');
                    legend('Either Strain','Non-AMR','AMR','fontsize',14,'Interpreter','latex');
                    box on;
                    grid on;
                    
             
                    
