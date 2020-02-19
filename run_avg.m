%% Simplified script to initialise, run and plot model output
% (single model run)

    clear all;
    clc;
    close all;
    
    NSims = 100;

    % General parmeters 
        N = 10000;          % population size
        n_Days = 10*365;     % days to simulate 
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
            params.p0 = [0.2 0.1 0];
        
    % display all parameters
        params
        
        
    
    %initialise counters    
    all_data = zeros(n_Days+1,2,10);
    all_vac_doses = zeros(n_Days+1,1);
    all_vac_current = zeros(n_Days+1,1);
    all_drug_doses = zeros(n_Days+1,2,1);
    %% initialise model (create new model object)
    for i = 1:NSims
        gono_model = VacAMR_IBM3(N, params, [], VERBOSE, LOW_MEM, [0,0,0]);
        %gono_model = VacAMR_IBM3(N, params, [], VERBOSE, LOW_MEM);
        
    %% run simulation for n_Days # of days
        %n_Days = 5*365;
        gono_model.simulate(n_Days);
        
    
    %% extract all counter data from model object 
    % (or can be referenced directly)
        data = gono_model.counters;
        prev_data = NSims*data.prevalence./N;
        all_data(:,:,i) = prev_data;
        all_vac_doses(:,i) = data.vac_doses_today;
        all_vac_current(:,i) = data.current_vac;
        all_drug_doses(:,1,i) = data.cipr;
        all_drug_doses(:,2,i) = data.cefta;
        
        %TODO avg drug doses
        %fprintf('Simulations complete = ')
        %reverseStr = '';
        %msg = [sprintf('Simulations complete... %3.2f', i)];
        %fprintf([reverseStr,msg]);
        %everseStr = repmat(sprintf('\b'), 1, length(msg))
        
        msg = ['Simulations completed = ',num2str(i), ' out of ', num2str(NSims)];
        disp(msg)
        
    end
        % plot prevalence time-series for whole simulation
        % (using my built in function)
            %gono_model.plot_prev(data, [0 n_Days], [])
        
        % manual plots (examples)
            % get prevalence (per strain) from counter varfiable
            % (not yet normalised with respect to the population size)
            
            
            %prev_data = 100*data.prevalence./N;
            %size(all_data)
            plot_data = mean(all_data,3);
            size(plot_data);
            %standard deviation at all points for plotting confidence
            %intervals, if we want to do that
            conf = std(all_data,0,3);
            size(conf)
            %quartiles for error bars that dont go below 0
            i25 = quantile(all_data,0.25,3); %25th percentile
            i75 = quantile(all_data,0.75,3);%75th percentile
            
            avg_vac_current = mean(all_vac_current,2);
            std_vac_current = std(all_vac_current,0,2);
            avg_vac_doses = mean(all_vac_doses,2);
            std_vac_doses = std(all_vac_doses,0,2);
            
            
           
           
            %plot 25th and 75th percentiles, not standard deviation as we
            %arent going to have negative prevalence so normal distribution
            %probably doesnt make sense here
            figure('name', 'Average Prevalence with interquartile range');
                hold on;
                %plot([0:n_Days],plot_data(:,1),'b-');
                %plot([0:n_Days],plot_data(:,2),'r-');
                shadedErrorBar([0:n_Days],plot_data(:,1),[conf(:,1),conf(:,1)],'lineprops','b');
                shadedErrorBar([0:n_Days],plot_data(:,2),[conf(:,2),conf(:,2)],'lineprops','r');
                %shadedErrorBar([0:n_Days],plot_data(:,1),[plot_data(:,1)-i25(:,1),i75(:,1)-plot_data(:,1)],'lineprops','b');
                %shadedErrorBar([0:n_Days],plot_data(:,2),[plot_data(:,2)-i25(:,2),i75(:,2)-plot_data(:,2)],'lineprops','r');
                legend('non-AMR','AMR');
                xlabel('Time (days)')
                ylabel('Average Prevalence (%)');
                box on;
                grid on;
            
            
            figure('name', 'Average Prevalence');
                hold on;
                plot([0:n_Days],plot_data(:,1),'b-');
                plot([0:n_Days],plot_data(:,2),'r-');
                legend('non-AMR','AMR');
                xlabel('Time (days)')
                ylabel('Average Prevalence (%)');
                box on;
                grid on;
            
            
                figure('name','Strain prevalence');
                    hold on;
                    plot([0:n_Days],prev_data(:,1),'b-'); % non AMR strain
                    plot([0:n_Days],prev_data(:,2),'r-'); % AMR strain
                    legend('non-AMR','AMR');
                    xlabel('Time (days)')
                    ylabel('Prevalence (%)');
                    box on;
                    grid on;
            
            % drug administration of each drug given by the cumulative sum
            % of the daily dosage of each drug
            %TODO avg drug doses
                figure('name','Dosage','color','w');
                    hold on;
                    plot([0:n_Days], cumsum(data.cipr),'b-');
                    plot([0:n_Days], cumsum(data.cefta),'r-');
                    legend('Cipr/A','Ceft/A','location','northwest');
                    xlabel('Time (days)');
                    ylabel('Number of doses')
                    title('Cumulative drug doses administered');
                    box on;
                    grid on;
                    
                figure('name','Dosage','color','w');
                    hold on;
                    data.vac_doses_today;
                    data.births;
                    plot([0:n_Days], cumsum(data.vac_doses_today),'b-');
                    plot([0:n_Days], cumsum(data.births),'r-');
                    legend('vaccine doses','births');
                    xlabel('Time (days)');
                    ylabel('Number of vaccine doses')
                    title('Cumulative vaccine doses administered');
                    box on;
                    grid on;    
                
                figure('name','Average doses of vaccine');
                    hold on;
                    %plot([0:n_Days], cumsum(avg_vac_doses));
                    shadedErrorBar([0:n_Days],cumsum(avg_vac_doses),[std_vac_doses(:,1),std_vac_doses(:,1)])
                    xlabel('Time (days)');
                    ylabel('No of vaccine doses given');
                    title('average vaccine doeses given');
                    box on;
                    grid on;
                    
                figure('name','average no. of vaccinated people');
                    hold on;
                    plot([0:n_Days], avg_vac_current);
                    shadedErrorBar([0:n_Days],avg_vac_current,[std_vac_current(:,1),std_vac_current(:,1)])
                    xlabel('Time (days)');
                    ylabel('average no of people protected by the vaccine');
                    title('average no of people protected');
                    box on;
                grid on;
             
                    
