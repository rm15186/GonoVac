%% Simplified script to initialise, run and plot model output
% (single model run)

    clear all;
    clc;
    close all;

    % General parmeters 
        N = 100;          % population size
        n_Days = 1*365;     % days to simulate
    
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
        
    %% initialise model (create new model object)
        gono_model = AMR_IBM(N, params, [], VERBOSE, LOW_MEM);
        
    %% run simulation for n_Days # of days
        n_Days = 365;
        gono_model.simulate(n_Days);
    
    %% extract all counter data from model object 
    % (or can be referenced directly)
            data = gono_model.counters

        % plot prevalence time-series for whole simulation
        % (using my built in function)
            gono_model.plot_prev(data, [0 n_Days], [])
        
        % manual plots (examples)
            % get prevalence (per strain) from counter variable
            % (not yet normalised with respect to the population size)
                prev_data = 100*data.prevalence./N
                
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
