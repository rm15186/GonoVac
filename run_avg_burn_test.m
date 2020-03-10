%% for BC Simplified script to initialise, run and plot model output 

    %clear all;
    %clc;
    %close all;
    %for bc shall we do 500 simulations, burn in of x days, start with 10/%
    %prevalence
    NSims = 100;

    % General parmeters 
        N = 10000;          % population size
        n_Days = 4000; %10*365;     % days to simulate about 5 years
        %fewer than 10 days and the error bars behave strangely
    
        VERBOSE = true;
        LOW_MEM = false;
        
    % load preset simulation parameters from external file
        load('base_params.mat','params');
     
    % you can edit (baseline) parameters by overwriting the preset values, e.g.
        params.P_SYMPTOMS = 0.5;
        params.LAB_DELAY_MEAN = 12;
        % initial strain prevalence:
            % p0(1) = overall initial prevalence of gonorrhea (0.1 = 10%)
            % p0(2) = proportion of positive cases with AMR component
            % p0(3) = proportion of coinfection given AMR
            %params.p0 = [0.2 0.1 0];% original initial prevalence
            params.p0 = [0.1 0.4 0]; % burn in takes less time think its just as good?
            %params.p0 = [0.04 0.5 0];
        
    % display all parameters
        params;
        
        
    %for working out averages
    %initialise counters, set size here for speed    
    all_data = zeros(n_Days+1,2,NSims); %prevalence of the individual strains
    all_either = zeros(n_Days+1,1,NSims); % prevalence of either strain
    all_vac_doses = zeros(n_Days+1,NSims); 
    all_vac_current = zeros(n_Days+1,NSims);
    all_cipr_doses = zeros(n_Days+1,NSims);
    all_cefta_doses = zeros(n_Days+1,NSims);
    all_burn_in_prev = zeros(3001,2,NSims);
    all_burn_in_prev_either = zeros(3001,1,NSims);
    %% initialise model (create new model object)
    tic
    for i = 1:NSims
        gono_model = VacAMR_IBM3(N, params, [], VERBOSE, LOW_MEM, [0,0,0]);
        %gono_model = VacAMR_IBM3(N, params, [], VERBOSE, LOW_MEM);
        
    %% run simulation for n_Days # of days
        %n_Days = 5*365;
        gono_model.simulate(n_Days);
        
    
    %% extract all counter data from model object unsupress this for BC?
    % (or can be referenced directly)
        data = gono_model.counters;
        prev_data = 100*data.prevalence./N; %prevalence of both strains
        prev_either_data = 100*data.prev_either/N; %prevalence of either strain
        all_either(:,i) = prev_either_data;
        all_data(:,:,i) = prev_data;
        
        all_vac_doses(:,i) = data.vac_doses_today;
        all_vac_current(:,i) = data.current_vac;
        all_cipr_doses(:,i) = data.cipr;
        all_cefta_doses(:,i) = data.cefta;
        all_burn_in_prev(:,:,i) = 100*data.burn_in_prevalence./N; %normalise
        all_burn_in_prev_either(:,i) = 100*data.burn_in_prevalence_either/N;
        
        %msg = ['Simulations completed = ',num2str(i), ' out of ', num2str(NSims)];
        %disp(msg);
        
    end
    t=toc;
    msg1 = ['time elapsed to run ' num2str(NSims) ' = ' num2str(t)];
    disp(msg1)
    
    
    
    %how much time does one simulation take, multiply that by 500 to see
    %how long a big ass run is going to take
    
        % plot prevalence time-series for whole simulation
        % (using my built in function)
            %gono_model.plot_prev(data, [0 n_Days], [])
        
        % manual plots (examples)
            % get prevalence (per strain) from counter varfiable
            % (not yet normalised with respect to the population size)
            
            %% give us the data
            prev_data = 100*data.prevalence./N;
            plot_data = mean(all_data,3)
            plot_either = mean(all_either,3)
            plot_burn_in_prev_either = mean(all_burn_in_prev_either,3) %2001?
            plot_burn_in_prev = mean(all_burn_in_prev,3)
            
            std_burn_in_prev_either = std(all_burn_in_prev,0,3)
            std_burn_in_prev = std(all_burn_in_prev,0,3)
            
            %confidence intervals
            %conf = std(all_data,0,3);
            %confe = std(all_either,0,3); %confidence interval for either strain
            
            %quartiles for error bars that dont go below 0
            i25 = quantile(all_data,0.25,3) %25th percentile
            i75 = quantile(all_data,0.75,3)%75th percentile
            
            i25e = quantile(all_either,0.25,3)
            i75e = quantile(all_either,0.75,3)
  
            
            avg_vac_current = mean(all_vac_current,2)
            std_vac_current = std(all_vac_current,0,2)
            avg_vac_doses = mean(all_vac_doses,2)
            std_vac_doses = std(all_vac_doses,0,2)
            avg_cipr_doses = mean(all_cipr_doses,2)
            std_cipr_doses = std(all_cipr_doses,0,2)
            avg_cefta_doses = mean(all_cefta_doses,2)
            std_cefta_doses = std(all_cefta_doses,0,2)


%% all the plots are irrelevant for bc well put this in another file
%             %plot 25th and 75th percentiles, not standard deviation as we
%             %arent going to have negative prevalence so normal distribution
%             %probably doesnt make sense here
%             figure('name', 'Average Prevalence with interquartile range');
%                 hold on;
%                 %plot([0:n_Days],plot_data(:,1),'b-'); %no error bars
%                 %plot([0:n_Days],plot_data(:,2),'r-');
%                 %one standard deviation 
%                 shadedErrorBar([0:n_Days],plot_data(:,1),[conf(:,1),conf(:,1)],'lineprops','b'); 
%                 shadedErrorBar([0:n_Days],plot_data(:,2),[conf(:,2),conf(:,2)],'lineprops','r');
%                 %plot([0:n_Days],plot_either,'k');
%                 shadedErrorBar([0:n_Days],plot_either,[confe,confe],'lineprops','k')
%                 %QUARTILES
%                 %shadedErrorBar([0:n_Days,either_strain,[either_strain-i25e,i75e-either_strain]],'lineprops','k');
%                 %shadedErrorBar([0:n_Days],plot_data(:,1),[plot_data(:,1)-i25(:,1),i75(:,1)-plot_data(:,1)],'lineprops','b');
%                 %shadedErrorBar([0:n_Days],plot_data(:,2),[plot_data(:,2)-i25(:,2),i75(:,2)-plot_data(:,2)],'lineprops','r');
%                 legend('non-AMR','AMR');
%                 xlabel('Time (days)')
%                 ylabel('Average Prevalence (%)');
%                 box on;
%                 grid on;
%             
%             %add the overall prevalence and make these thicker
%             figure('name', 'Average Prevalence');
%                 hold on;
%                 plot([0:n_Days],plot_data(:,1),'b-');
%                 plot([0:n_Days],plot_data(:,2),'r-');
%                 plot([0:n_Days],plot_either,'k-');
%                 legend('non-AMR','AMR','Total Prevalence');
%                 xlabel('Time (days)')
%                 ylabel('Average Prevalence (%)');
%                 box on;
%                 grid on;
%             
%             
% %                 figure('name','Strain prevalence');
% %                     hold on;
% %                     plot([0:n_Days],prev_data(:,1),'b-'); % non AMR strain
% %                     plot([0:n_Days],prev_data(:,2),'r-'); % AMR strain
% %                     legend('non-AMR','AMR');
% %                     xlabel('Time (days)')
% %                     ylabel('Prevalence (%)');
% %                     box on;
% %                     grid on;
%                     
%                     
%             
%             % drug administration of each drug given by the cumulative sum
%             % of the daily dosage of each drug
%                 figure('name','Dosage','color','w');
%                     hold on;
%                     %plot([0:n_Days], cumsum(avg_cipr_doses));
%                     %plot([0:n_Days], cumsum(avg_cefta_doses));
% %                   size([0:n_Days]);
% %                   size(cumsum(avg_cipr_doses(:,1)));
% %                   size(std_cipr_doses(:,:,1)),
%                     shadedErrorBar([0:n_Days], cumsum(avg_cipr_doses(:,1)),[std_cipr_doses(:,:,1),std_cipr_doses(:,:,1)]);
%                     shadedErrorBar([0:n_Days],cumsum(avg_cefta_doses(:,1)),[std_cefta_doses(:,:,1),std_cefta_doses(:,:,1)]);
%                     legend('Cipr/A','Ceft/A','location','northwest');
%                     xlabel('Time (days)');
%                     ylabel('Number of doses')
%                     title('Cumulative drug doses administered');
%                     box on;
%                     grid on;
%                     
% %                 figure('name','Dosage','color','w');
% %                     hold on;
% %                     data.vac_doses_today;
% %                     data.births;
% %                     plot([0:n_Days], cumsum(data.vac_doses_today),'b-');
% %                     plot([0:n_Days], cumsum(data.births),'r-');
% %                     legend('vaccine doses','births');
% %                     xlabel('Time (days)');
% %                     ylabel('Number of vaccine doses')
% %                     title('Cumulative vaccine doses administered');
% %                     box on;
% %                     grid on;    
%                 
%                 figure('name','Average doses of vaccine');
%                     hold on;
%                     %plot([0:n_Days], cumsum(avg_vac_doses));
%                     shadedErrorBar([0:n_Days],cumsum(avg_vac_doses),[std_vac_doses(:,1),std_vac_doses(:,1)])
%                     xlabel('Time (days)');
%                     ylabel('No of vaccine doses given');
%                     title('average vaccine doeses given');
%                     box on;
%                     grid on;
%                     
%                 figure('name','average no. of vaccinated people');
%                     hold on;
%                     plot([0:n_Days], avg_vac_current);
%                     shadedErrorBar([0:n_Days],avg_vac_current,[std_vac_current(:,1),std_vac_current(:,1)])
%                     xlabel('Time (days)');
%                     ylabel('average no of people protected by the vaccine');
%                     title('average no of people protected');
%                     box on;
%                 grid on;
%                 
%                 %TODO plot end of burn in and then plot effect of vaccine
%                 %currently not working but should just plot what happens
%                 %over the 2000 day burn in 
%                  figure('name','burn in');
%                      hold on;
%                      plot([0:3000],plot_burn_in_prev_either,'k');
%                      plot([0:3000],plot_burn_in_prev(:,1),'b');
%                      plot([0:3000],plot_burn_in_prev(:,2),'r');
%                      xlabel('Time (days)');
%                      ylabel('prevalence in burn in');
%                      title('prevalence over burn in period');
%                      box on;
%                      grid on;
%                      
%                  
%                  %stitch the data together 
%                  days = [0:500+n_Days]; %end of burn in then length of simulation
%                  burn = plot_burn_in_prev_either([2500:2999]);
%                  %plot_change = cat(1,burn,plot_data);
%                  plot_change = cat(1,burn,plot_either);
%                  size(plot_change)
%                  size(days)
%                  figure('name','steady state then vaccine');
%                     hold on;
%                     plot(days,plot_change);
%                     xline(500);
%                     xlabel('Time (days)');
%                     ylabel('Prevalence (%)');
%                     title('Impact of vaccine on prevalence');
%                     box on;
%                     grid on;
%                     
%              
%                     
