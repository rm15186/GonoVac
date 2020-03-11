%% prevalence
plot_data = mean(all_data,3);
            plot_either = mean(all_either,3);
            
            
            %standard deviation at all points for plotting confidence
            %intervals, if we want to do that
            conf = std(all_data,0,3);
            confe = std(all_either,0,2);
            %size(conf);
            %quartiles for error bars that dont go below 0
            i25 = quantile(all_data,0.25,3); %25th percentile
            i75 = quantile(all_data,0.75,3);%75th percentile
            
            i25e = quantile(all_either,0.25,3);
            i75e = quantile(all_either,0.75,3);
  
            %%burn in prevalence 
            plot_burn_in_prev_either = mean(all_burn_in_prev_either,3); %2001?
            plot_burn_in_prev = mean(all_burn_in_prev,2);
            std_burn_in_prev_either = std(all_burn_in_prev_either,0,3);
            std_burn_in_prev = std(all_burn_in_prev,0,2);
            
            %% current people vaccinated %
            avg_vac_current = mean(all_vac_current,2);
            std_vac_current = std(all_vac_current,0,2);
            range_vac_current = [quantile(all_vac_current,0.25,2),quantile(all_vac_current,0.75,2)];
            plot_range_vac_current = [avg_vac_current-range_vac_current(:,1),range_vac_current(:,2)-avg_vac_current];
            %% doses of vaccine given
            avg_vac_doses = mean(all_vac_doses,2);
            std_vac_doses = std(all_vac_doses,0,2);
            range_vac_doses = [quantile(all_vac_doses,0.25,2),quantile(all_vac_doses,0.75,2)];
            plot_range_vac_doses = [avg_vac_doses-range_vac_doses(:,1),avg_vac_doses-range_vac_doses(:,2)-avg_vac_doses];
            
            %% doses of antibiotics - cipr not used 
            avg_cipr_doses = mean(all_cipr_doses,2);
            std_cipr_doses = std(all_cipr_doses,0,2);
            
            avg_cefta_doses = mean(all_cefta_doses,2);
            std_cefta_doses = std(all_cefta_doses,0,2);
            range_cefta_doses = [quantile(all_cefta_doses,0.25,2),quantile(all_cefta_doses,0.75,2)];
            plot_range_cefta_doses = [avg_cefta_doses-range_cefta_doses(:,1),range_cefta_doses(:,2)-avg_cefta_doses];

            
            size(all_cipr_doses);
            size(mean(all_cipr_doses,2));
            
            %% plot
            close all
figure('name', 'Average Prevalence with interquartile range');
                hold on;
                
                size(0:n_Days);
                size(plot_either);
                size([plot_either-i25e,i75e-plot_either]);
                shadedErrorBar([0:n_Days],plot_either,[plot_either-i25e,i75e-plot_either],'lineprops','k');
                shadedErrorBar([0:n_Days],plot_data(:,1),[plot_data(:,1)-i25(:,1),i75(:,1)-plot_data(:,1)],'lineprops','b');
                shadedErrorBar([0:n_Days],plot_data(:,2),[plot_data(:,2)-i25(:,2),i75(:,2)-plot_data(:,2)],'lineprops','r');
                legend('Either Strain','non-AMR','AMR');
                xlabel('Time (days)')
                ylabel('Average Prevalence (%)');
                box on;
                grid on;
            
            %add the overall prevalence and make these thicker
            figure('name', 'Average Prevalence');
                hold on;
                plot([0:n_Days],plot_data(:,1),'b-');
                plot([0:n_Days],plot_data(:,2),'r-');
                plot([0:n_Days],plot_either,'k-');
                legend('non-AMR','AMR','Total Prevalence');
                xlabel('Time (days)')
                ylabel('Average Prevalence (%)');
                box on;
                grid on;
                         
            
            % drug administration of each drug given by the cumulative sum
            % of the daily dosage of each drug
                figure('name','Dosage','color','w');
                    hold on;
                    %shadedErrorBar([0:n_Days], cumsum(avg_cipr_doses(:,1)),[std_cipr_doses(:,:,1),std_cipr_doses(:,:,1)]);
                    shadedErrorBar([0:n_Days],cumsum(avg_cefta_doses(:,1)),plot_range_cefta_doses);
                    %legend('Cipr/A','Ceft/A','location','northwest');
                    xlabel('Time (days)');
                    ylabel('Number of doses')
                    title('Cumulative drug doses administered');
                    box on;
                    grid on;
                    

                
                figure('name','Average doses of vaccine');
                    hold on;
                    %plot([0:n_Days], cumsum(avg_vac_doses));
                    shadedErrorBar([0:n_Days],cumsum(avg_vac_doses),plot_range_vac_doses)
                    xlabel('Time (days)');
                    ylabel('No of vaccine doses given');
                    title('average vaccine doeses given');
                    box on;
                    grid on;
                    
                figure('name','average no. of vaccinated people');
                    hold on;
                    plot([0:n_Days], avg_vac_current);
                    shadedErrorBar([0:n_Days],avg_vac_current,plot_range_vac_current)
                    xlabel('Time (days)');
                    ylabel('average no of people protected by the vaccine');
                    title('average no of people protected');
                    box on;
                grid on;
                
                %TODO plot end of burn in and then plot effect of vaccine
                %currently not working but should just plot what happens
                %over the 2000 day burn in 
                 figure('name','burn in');
                     hold on;
                     plot([0:3000],plot_burn_in_prev_either,'k');
                     plot([0:3000],plot_burn_in_prev(:,1),'b');
                     plot([0:3000],plot_burn_in_prev(:,2),'r');
                     shadedErrorBar([0:3000],plot_burn_in_prev_either,[std_burn_in_prev_either],'lineprops','k')
                     shadedErrorBar([0:3000],plot_burn_in_prev(:,1),[std_burn_in_prev(:,1),std_burn_in_prev(:,1)],'lineprops','b')
                     shadedErrorBar([0:3000],plot_burn_in_prev(:,2),[std_burn_in_prev(:,2),std_burn_in_prev(:,2)],'lineprops','r')
                     xlabel('Time (days)');
                     ylabel('prevalence in burn in');
                     title('prevalence over burn in period');
                     box on;
                     grid on;
                     
                 
                 %stitch the data together 
                 days = [0:500+n_Days]; %end of burn in then length of simulation
                 burn = plot_burn_in_prev_either([2500:2999]);
                 burn_strains = plot_burn_in_prev([2500:2999],:);
                 %plot_change = cat(1,burn,plot_data);
                 plot_change = cat(1,burn,plot_either);
                 plot_change_both = cat(1,burn_strains,plot_data);
                 size(plot_change);
                 size(days);
                 figure('name','steady state then vaccine');
                    hold on;
                    plot(days,plot_change,'k');
                    plot(days,plot_change_both(:,1),'b');
                    plot(days,plot_change_both(:,2),'r');
                    xline(500);
                    xlabel('Time (days)','Interpreter','latex');
                    ylabel('Prevalence (%)','Interpreter','latex');
                    legend('Either Strain','Non-AMR','AMR','fontsize',14,'Interpreter','latex');
                    title('Impact of vaccine on prevalence','Interpreter','latex');
                    box on;
                    grid on;