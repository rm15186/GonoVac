      
                figure('name','average no. of vaccinated people');
                    hold on;
                   % avg_vac_current3High = 100*avg_vac_current/N;
                    %std_vac_current = 100^2*std_vac_current/N^2; %normalise
                    
                    plot([0:n_Days],vac_people31High,'linewidth',2);
                    %shadedErrorBar([0:n_Days], vac_people31, [plot_range31High(:,1),plot_range31High(:,2)])
                    plot([0:n_Days], vac_people62High,'linewidth',2);
                    %shadedErrorBar([0:n_Days], vac_people31, [plot_range_vac_current(:,1),plot_range_vac_current(:,2)])
                    plot([0:n_Days], vac_people100High,'linewidth',2);
                    %shadedErrorBar([0:n_Days], vac_people31, [plot_range_vac_current(:,1),plot_range_vac_current(:,2)])
                    
                    %shadedErrorBar([0:n_Days],avg_vac_current,[std_vac_current(:,1),std_vac_current(:,1)])
                    xlabel('Time (days)','fontsize',22,'Interpreter','latex');
                    ylabel('People protected (\%)','fontsize',22,'Interpreter','latex');
                    legend('31\%','62\%','100\%','fontsize',20,'Interpreter','latex');
                    %title('average no of people protected','fontsize',14,'Interpreter','latex');
                    box on;
                    grid on;
                    
                    
%         figure('name','average no. of vaccinated people');
%             hold on;
%            % avg_vac_current3High = 100*avg_vac_current/N;
%             %std_vac_current = 100^2*std_vac_current/N^2; %normalise
% 
%             plot([0:n_Days], vac_people31,'linewidth',2);
%             %shadedErrorBar([0:n_Days], vac_people31, [plot_range_vac_current(:,1),plot_range_vac_current(:,2)])
%             plot([0:n_Days], vac_people62,'linewidth',2);
%             %shadedErrorBar([0:n_Days], vac_people31, [plot_range_vac_current(:,1),plot_range_vac_current(:,2)])
%             plot([0:n_Days], vac_people100,'linewidth',2);
%             %shadedErrorBar([0:n_Days], vac_people31, [plot_range_vac_current(:,1),plot_range_vac_current(:,2)])
% 
%             %shadedErrorBar([0:n_Days],avg_vac_current,[std_vac_current(:,1),std_vac_current(:,1)])
%             xlabel('Time (days)','fontsize',22,'Interpreter','latex');
%             ylabel('People protected (\%)','fontsize',22,'Interpreter','latex');
%             legend('31\%','62\%','100\%','fontsize',22,'Interpreter','latex');
%             %title('average no of people protected','fontsize',14,'Interpreter','latex');
%             box on;
%             grid on;