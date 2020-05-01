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
        gono_model = VacAMR_IBM3(N, params, [], VERBOSE, LOW_MEM, [0,1,0]);
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
    msg1 = ['time elapsed to run ' num2str(NSims) ' simulations = ' num2str(t)];
    disp(msg1)
    msg2 = '];';
    
            plot_data = mean(all_data,3)
            disp(msg2)
            plot_either = mean(all_either,3) 
            disp(msg2)
            
            %standard deviation at all points for plotting confidence
            %intervals, if we want to do that
            conf = std(all_data,0,3)
            disp(msg2)
            confe = std(all_either,0,3)
            disp(msg2)
            %size(conf);
            %quartiles for error bars that dont go below 0
            i25 = quantile(all_data,0.25,3) %25th percentile
            disp(msg2)
            i75 = quantile(all_data,0.75,3)%75th percentile
            disp(msg2)
            
            i25e = quantile(all_either,0.25,3)
            disp(msg2)
            i75e = quantile(all_either,0.75,3)
            disp(msg2)
            
            %%burn in prevalence 
            plot_burn_in_prev_either = mean(all_burn_in_prev_either,3) %2001? %10
            disp(msg2) 
            plot_burn_in_prev = mean(all_burn_in_prev,3) %6,4
            disp(msg2)
            std_burn_in_prev_either = std(all_burn_in_prev_either,0,3)
            disp(msg2)
            std_burn_in_prev = std(all_burn_in_prev,0,3)
            disp(msg2)
            
            %TODO write something that plots this! NEW delete this if it
            %breaks on 50 not 45
            %these are going negative and i might jjust not bother
            i25burn1 = quantile(all_burn_in_prev(:,1),0.25,3)
            disp(msg2)
            i25burn2 = quantile(all_burn_in_prev(:,2),0.25,3)
            disp(msg2)
            i75burn1 = quantile(all_burn_in_prev(:,1),0.75,3)
            disp(msg2)
            i75burn2 = quantile(all_burn_in_prev(:,2),0.75,3)
            disp(msg2)
            i25burn_either = quantile(all_burn_in_prev_either,0.25,3)
            disp(msg2)
            i75burn_either = quantile(all_burn_in_prev_either,0.75,3)
            disp(msg2)
            plot_range_burn_either=[plot_burn_in_prev_either-i25burn_either,i75burn_either-plot_burn_in_prev_either]
            disp(msg2)
            plot_burn_1 = [plot_burn_in_prev(:,1)-i25burn1,i75burn1-plot_burn_in_prev(:,1)]
            disp(msg2)
            plot_burn_2 = [plot_burn_in_prev(:,2)-i25burn2,i75burn2-plot_burn_in_prev(:,2)]
            disp(msg2)
            

            %% current people vaccinated %
            avg_vac_current = mean(all_vac_current,2)
            disp(msg2)
            std_vac_current = std(all_vac_current,0,2)
            disp(msg2)
            range_vac_current = [quantile(all_vac_current,0.25,2),quantile(all_vac_current,0.75,2)]
            disp(msg2)
            plot_range_vac_current = [avg_vac_current-range_vac_current(:,1),range_vac_current(:,2)-avg_vac_current]
            disp(msg2)
            %% doses of vaccine given
            avg_vac_doses = mean(all_vac_doses,2)
            disp(msg2)
            std_vac_doses = std(all_vac_doses,0,2)
            disp(msg2)
            range_vac_doses1 = quantile(all_vac_doses,0.25,2)
            disp(msg2)
            range_vac_doses2 = quantile(all_vac_doses,0.75,2)
            disp(msg2)
            plot_range_vac_doses = [avg_vac_doses-range_vac_doses1,range_vac_doses2-avg_vac_doses]
            disp(msg2)
            %% doses of antibiotics - cipr not used 
            avg_cipr_doses = mean(all_cipr_doses,2)
            disp(msg2)
            std_cipr_doses = std(all_cipr_doses,0,2)
            disp(msg2)
            
            avg_cefta_doses = mean(all_cefta_doses,2)
            disp(msg2)
            std_cefta_doses = std(all_cefta_doses,0,2)
            disp(msg2)
            range_cefta_doses = [quantile(all_cefta_doses,0.25,2),quantile(all_cefta_doses,0.75,2)]
            disp(msg2)
            plot_range_cefta_doses = [avg_cefta_doses-range_cefta_doses(:,1),range_cefta_doses(:,2)-avg_cefta_doses]
            disp(msg2)
            
            a = zeros(NSims,1);
            for i = 1:NSims
                a(i) = all_either(4000,1,i);
            end
            finalPrev = a
            disp(']')



    
    