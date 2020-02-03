%% RM adding vaccination into model via a new vector
%% Two-strain (nonAMR/AMR) coinfection SIS model for Gonorrhea
%  This is the main class library - run model via the run_simple.m script
%   
%   - Static, scale-free partnership network
%
%   - Independent transmission of strains with rate BETA_[1,2]
%
%   - Natural recovery with rate R (both strains independently)
%
%   - Recovery via birth/death     (all infections removed)
%
%   - Recovery via intervention:   - voluntary treatment seeking 
%                                   : for symptomatic patients
%                                    (a proportion of which seek care)
%
%                                  - screening with rate GAMMA
%                                   : updates RECALL notifications for infectees
%
%                                  - partner tracing 
%                                   : proportion PSI of infectees partners
%                                   are traced with the correct nonAMR/AMR
%                                   status

%
%
%
%  AUTHOR:  Adam Zienkiewicz, 2016 (adamziek@gmail.com)
%           University of Bristol
%
%  LAST UPDATED:   20/10/16

classdef AMR_IBM < handle

    properties
        
        % Sexual partnership network parameters
        % --------------------------------      
        N;              % population size (number of individuals)   
        ALPHA;          % power-law network (slope) P(k) ~ k^-ALPHA
        FULL_MAX_PARTNERS; % max. number of partners
        
        % SIS parameters (same for all strains)
        % ----------------------------------------
        n_Strains;      % number of different disease strains 
        R;              % natural recovery rate (/day)
        MU;             % birth / death rate (/day)
        BETA;           % transmission rate (/day)
        GAMMA;          % screening rate (/day)
        PSI;            % tracing efficiency
        MAX_TRACE;      % max. number of partners who can be traced 
                        %   per infected index patient
        
        % initial infection prevalence
        % p0(1) - overall prevalence (either strain)
        % p0(2) - proportion of AMR to nonAMR strain
        %   i.e.  
        %        if 40% of all gonorrhea cases are found to be AMR 
        %        then p0(2) = 0.4
        p0;             
        
        P_SYMPTOMS;          % proportion of infected cases which are symptomatic
        P_SEEKS_TREATMENT;   % proportion of symptomatic individuals who will seek treatment
        
        % maximum delay between begin flagged as nonAMR (recall/trace) and being treated as nonAMR
        NON_AMR_MAX_DELAY;
        
        % maximum # partners for a recalled/traced individual to be treated as nonAMR 
        NON_AMR_MAX_PARTNERS;
        
        RESTRICT_MAX_PARTNERS;  % max. # partners in a given time-window
        RESTRICT_RATE;          % time-window duration for restricted partnership network (days)
        
        % Time delays (indays) for
        % lab results / symptom onset / seek / recall and traces
        LAB_DELAY_MEAN;     % mean delay between being screened/treated and AMR risk being known
        LAB_DELAY_STD;      % std. of above
        ONSET_DELAY_MEAN;   % mean delay between becoming infected and showing symptoms (if any)
        ONSET_DELAY_STD;    % std. of above
        SEEK_DELAY_MEAN;    % mean number of days between symptom onset and treatment seeking (arrival)
        SEEK_DELAY_STD;     % std. of above for normal dist.
        RECALL_DELAY_MEAN;  % mean number of days between recall notication and arrival
        RECALL_DELAY_STD;   % std. of above for normal dist.
        TRACE_DELAY_MEAN;   % mean number of days between trace notication and arrival
        TRACE_DELAY_STD;    % std. of above for normal dist.   
        
        % proportion currently treated with recommended dual4
        % therapy (Ceft/A) at point of care, without prior
        % knowledge of strain susceptibility
        P_BLINDTREAT_AS_AMR;
        
        % allow the extra nonAMR treatment route using recall/trace data
        % true / false
        % (nonAMR treatment can still occur if P_BLINDTREAT_AS_AMR < 1)
        ENABLE_nonAMR_RECALL;
        ENABLE_nonAMR_TRACE;
        
        % require all traced individuals to be screened rather than treated
        % on appointment (true/false)
        PRESCREEN_TRACED;
        
        % require all symptomatic voluntary treatment seeking individuals 
        % to be screened rather than treated on appointment (true/false)
        PRESCREEN_SEEKED;
        
        % Enable POCT test for all treatment seekers (all routes)
        ENABLE_POCT;
        
        % Is POCT strain discriminatory - or 100% Ceft (FALSE)
        DISCRIM_POCT;
        
        % if ALLOW_COINFECTION = true then individuals can be infected with
        % both nonAMR and AMR strains simulteneously.
        % if false, then an individual already infected with one strain,
        % cannot contract the other
        ALLOW_COINFECTION;
        
        % ENABLE / DISABLE any form of treatment
        ALLOW_TREAT;
       

        % network toplogy
        % ---------------
        adj;                    % current sexual partnership network (adjacency matrix)
        adj_full;               % partnership network including all relationships over a year
        adj_int;                % integrated partnership network from day zero to today
        rel_list_full;          % undirected list of relationships (node-node))
        num_partners_full;      % degree sequence of full (yearly) network
        mean_degree_full;       % mean degree (ave.# partners)
        id_casuals;             % list of which nodes will require pruning / number of links to prung
        fixed_rels;
        other_rels;
        
        adj_xmin;
        adj_L;
        
        num_partners;           % degree sequence of current partnership network
        mean_degree_current;    % mean degree of current partnership network (stored history)
        n_Comp;                 % number of connected network components
        comp_sizes;             % sizes of connected components
        comp_members;           % members of connected components
        
        % data arrays and counters (infection state etc.)
        % ----------------------------------------------------
        indiv_state;        % state matrix (individuals infected by each strain over time)
        infected_since;     % day last infection acquired (per stra
        infection_count;    % number of times individual has had infection (per strain)
        symptoms;           % whether current infection is symptomatic (single value)
        seeks_treatment_on; % day a symptomatic patient will seek treatment (without waiting for recall/trace)
        
        % screening (recall) and tracing notification arrays
        % Nx3 arrays 
        %  col 1: day of recall/trace issue
        %  col 2: treatement day
        %  col 3: AMR risk flag
        recall_notify;
        trace_notify;
        
        today;              % number of days which have been simulated
        
        burn_in;            % whether or not we are in a 'burn-in' cycle
                            %   to achieve a set AMR ratio etc.
        
        % counters
        counters;           % structure array of counters (data outputs)
%       with fields:        
%             counters.drug_count;         % # times each nonAMR/AMR drug prescribed to each patient
%             counters.cefta;              % # CEFT/A (AMR) prescriptions on each day of simulation
%             counters.cefta_notinf;       % # CEFT/A prescriptions which weren't required (susceptible)
%             counters.cefta_nonAMR;       % # CEFT/A prescriptions which weren't required (nonAMR)
%             counters.cefta_AMR;          % # CEFT/A correct prescriptions (AMR infected)
%             counters.cipr;               % # Cipr (nonAMR) prescriptions on each day of simulation
%             counters.cipr_notinf;        % # Cipr prescriptions which weren't required (susceptible)
%             counters.cipr_nonAMR;        % # Cipr correct prescriptions (only nonAMR infected)
%             counters.cipr_AMR;         % # Cipr prescriptions which were incorrect (leading to recalled AMR patient)
%             counters.prevalence;         % prevalences - number of individuals with each strain over time
%             counters.prev_both;          % prevalences - number of individuals with both strains (coinfected) over time
%             counters.prev_either;    int8    % prevalence (either strain or both)
%             counters.incidence;          % infection incidence (number of new infections each day, per strain)
%             counters.incidence_either;   % infection incidence (new infections of either strain each day)
%             counters.n_attended_seeked;  % number of attendees who voluntarily seeked treatment each time step
%             counters.n_attended_recalled;% number of attendees who were recalled each time step
%             counters.n_attended_traced;  % number of attendees who were traced each time step
%             counters.n_attended_treated; % number of attendees who were actually treated each time step
%             counters.n_screened;         % number of people screened at each time step
%             counters.n_traced;           % number of people traced at each time step

        ipr;    % infected partner ratio (per node degree)

        % Other
        % ------
        param_updates;  % list of commands and dates to update parameters during simulation 
        DIR_NAME;  % directory name for saving data files and plots
        fig_h;     % figure handles
        plot_names; % name of plot (for saviandng purposes)
        
        VERBOSE; % for controlling console text output
        LOW_MEM; % enable low memory mode (no state histories of each individual)
        DEBUG;   % extra console output for debugging
        
    end % class properties

        methods
            % Class Constructor
            function [self] = AMR_IBM(N, params, adj_set, VERBOSE, LOW_MEM)
                
                % switch on/off console output if specified
                if nargin == 4
                    self.VERBOSE = VERBOSE;
                elseif nargin == 5
                    self.VERBOSE = VERBOSE;
                    self.LOW_MEM = LOW_MEM;
                end
                
                % debug console output ON/OFF
                self.DEBUG = false;
                
                %  set values from inputs
                self.N = N;
                
                % for this AMR / non-AMR implementation, fix n_Strains=2
                % (dont even attempt changing this value as much of the code now
                % relies on this being equal to 2!)
                self.n_Strains = 2;
                
                % Set parameter values from params structure
                if isempty(params)
                  error('No model parameters have been assigned!');
                end   
                self.ALPHA = params.ALPHA;
                self.FULL_MAX_PARTNERS = params.FULL_MAX_PARTNERS;
                self.RESTRICT_MAX_PARTNERS = params.RESTRICT_MAX_PARTNERS;
                self.RESTRICT_RATE = params.RESTRICT_RATE;
                self.R = params.R;
                self.MU = params.MU;
                self.BETA = params.BETA;
                self.GAMMA = params.GAMMA;
                self.PSI = params.PSI;
                self.MAX_TRACE = params.MAX_TRACE;
                self.LAB_DELAY_MEAN = params.LAB_DELAY_MEAN;
                self.LAB_DELAY_STD = params.LAB_DELAY_STD;
                self.ONSET_DELAY_MEAN = params.ONSET_DELAY_MEAN;
                self.ONSET_DELAY_STD = params.ONSET_DELAY_STD;
                self.SEEK_DELAY_MEAN = params.SEEK_DELAY_MEAN;
                self.SEEK_DELAY_STD = params.SEEK_DELAY_STD;
                self.RECALL_DELAY_MEAN = params.RECALL_DELAY_MEAN;
                self.RECALL_DELAY_STD = params.RECALL_DELAY_STD;
                self.TRACE_DELAY_MEAN = params.TRACE_DELAY_MEAN;
                self.TRACE_DELAY_STD = params.TRACE_DELAY_STD;
                self.P_SYMPTOMS = params.P_SYMPTOMS;  
                self.P_SEEKS_TREATMENT = params.P_SEEKS_TREATMENT;
                self.NON_AMR_MAX_DELAY = params.NON_AMR_MAX_DELAY;
                self.NON_AMR_MAX_PARTNERS = params.NON_AMR_MAX_PARTNERS;
                self.P_BLINDTREAT_AS_AMR = params.P_BLINDTREAT_AS_AMR; 
                self.ENABLE_nonAMR_RECALL = params.ENABLE_nonAMR_RECALL;
                self.ENABLE_nonAMR_TRACE = params.ENABLE_nonAMR_TRACE;
                self.ENABLE_POCT = params.ENABLE_POCT;
                self.DISCRIM_POCT = params.DISCRIM_POCT;
                self.PRESCREEN_TRACED = params.PRESCREEN_TRACED;
                self.PRESCREEN_SEEKED = params.PRESCREEN_SEEKED;
                self.ALLOW_COINFECTION = params.ALLOW_COINFECTION;
                self.ALLOW_TREAT = params.ALLOW_TREAT;
       
                %-----------------------------------------------------------------
                
                
                % some warnings / overrides
                if self.ENABLE_POCT && ((self.LAB_DELAY_MEAN > 0) || (self.LAB_DELAY_STD > 0))
                    warning('Point-of-care treatment has been enabled but LAB delay values are non-zero...resetting delays')
                    self.LAB_DELAY_MEAN = 0;
                    self.LAB_DELAY_STD = 0;
                end
  
                % set initial prevalences
                self.p0 = params.p0;
                               
                % state matrix 
                % 1: non-AMR state
                % 2: AMR state
                % 3: date infection acquired (either strain)
                % 4: symptomatic (true / false)
                self.indiv_state = false(self.N,self.n_Strains,1); 
                size(self.indiv_state)
                
                % day on which particular strain passed to (susceptible) individual
                % initialise with NaNs
                self.infected_since = nan(self.N,self.n_Strains);
                    
                % flags indicating whether infections are symptomatic
                self.symptoms = false(self.N,1);
                
                % values indicating if/when(day) symptomatic individual will voluntarily seek
                % treatment (init with NaNs - value assigned when necessary)
                self.seeks_treatment_on = nan(self.N,1);
                
                % screening / recall notification matrix
                % (day recall issued, day attending, AMR flag)
                self.recall_notify = nan(self.N,3);
                
                % tracing notification matr2000ix
                % (day trace issued, day attending, AMR flag)
                self.trace_notify = nan(self.N,3);
                
                % treatment counters for each individual
                self.counters = struct;
                
                    self.counters.N = self.N;
                    
                    % count number of times each individual has been infected
                    % (per strain)
                    self.counters.infection_count = zeros(self.N,self.n_Strains);
                
                    % 2 drugs,e.g. Cipr for non-AMR, Ceft/A for AMR strains
                    self.counters.drug_count = zeros(self.N,2); 

                    % total CEFT/A prescribed each day
                    self.counters.cefta = 0; 
                    
                    % of which prescribed to non-infected individuals [0 0]
                    self.counters.cefta_notinf = 0;
                    
                    % of which prescribed to nonAMR individuals [1 0]
                    % (POOR TREATMENT CHOICE / WASTE)                            
                    self.counters.cefta_nonAMR = 0;
                    
                    % of which prescribed to AMR individuals [0 1] or [1 1]
                    % (CORRECT TREATMENT)
                    self.counters.cefta_AMR = 0;
                    
                    % total Cipr prescribed each day
                    self.counters.cipr = 0; 
                    
                    % of which prescribed to non-infected individuals [0 0]
                    self.counters.cipr_notinf = 0;
                    
                    % of which prescribed to AMR infected individuals [0 1] or [1 1 ]
                    % (MISTREATMENT REQUIRING RECALL)
                    self.counters.cipr_AMR = 0; % each day
                    
                    % of which prescribed to AMR infected individuals [1 0]
                    % (CORRECT TREATMENT)
                    self.counters.cipr_nonAMR = 0; % each day
                    
                    % prevalences - number of individuals with each strain
                    % : nonAMR / AMR
                    self.counters.prevalence = zeros(1,2);
                    self.counters.prev_both = 0;
                    self.counters.prev_either = 0;
                    self.counters.prev_nonAMRonly = 0;

                    % incidence - number of new infections of each strain each day
                    self.counters.incidence = zeros(1,2);
                    self.counters.incidence_either = 0;
                    self.counters.incidence_new = 0;
                    
                    %TODO RM add vaccinated counter
                    

                    % number of people who attended after voluntarily seeking at each time step
                    % (2nd/3rd column: number infected with each strain)
                    self.counters.n_attended_seeked = zeros(1,3);
                    
                    % number of people who attended after being recalled at each time step
                    % (2nd/3rd column: number infected with each strain)
                    self.counters.n_attended_recalled = zeros(1,3);
                    
                    % number of people who attended after being recalled at each time step
                    % (2nd/3rd column: number infected with each strain)
                    self.counters.n_attended_traced = zeros(1,3);
                    
                    % number of people who attended and were treated at each time step
                    % (2nd/3rd column: number infected with each strain)
                    self.counters.n_attended_treated = zeros(1,3);
                    
                    self.counters.n_treated_infected_symptomatic = zeros(1,2);
                    
                    % number of people who were screened at each time step
                    % (2nd/3rd column: number infected with each strain)
                    % - number of recalled individuals will be rowsum of
                    % columns 2 and 3
                    self.counters.n_screened = zeros(1,3);

                    % number of people traced at each time step
                    % (2nd/3rd column: number infected with each strain)
                    self.counters.n_traced = zeros(1,3);
                    
                    % number of people recalled at each time step
                    % (2nd/3rd column: number infected with each strain)
                    self.counters.n_recalled = zeros(1,3);
                    
                    % data for analysing undertreatment with nonAMR
                    % therapy using informed treatment scenario
                    self.counters.informed.undertreated_degree = [];
                    self.counters.informed.undertreated_infected_duration = [];
                    self.counters.informed.correct_degree = [];
                    self.counters.informed.correct_infected_duration = [];
                
                % file info for saving 
                %self.DIR_NAME = ['_N=',num2str(N),'_nS=',num2str(self.n_Strains),'_days=',num2str(self.today)];
                self.DIR_NAME = 'temp';
                self.fig_h = [];
                self.plot_names = {};
                self.param_updates = {};
                    
                % console output
                if self.VERBOSE
                    fprintf(['\n--------------------------------------------------\n',...
                                'Multi-strain SIS\n',...
                                '------------------\n']);
                    fprintf(['Population size (N) \t\t= ',num2str(self.N),...
                        '\nNumber of Strains \t\t= ',num2str(self.n_Strains),...
                        '\nInfection rate (Beta) \t\t= ',sprintf('%.2s',self.BETA(1)),' (non-AMR), ',sprintf('%.2s',self.BETA(2)),' (AMR)',...
                        '\nNatural recovery rate (R) \t= ',sprintf('%.2s',self.R),...
                        '\nBirth/death rate (MU) \t\t= ',sprintf('%.2s',self.MU),...
                        '\nScreening rate (GAMMA) \t\t= ',sprintf('%.2s',self.GAMMA),...
                        '\nTracing efficiency (PSI) \t= ',sprintf('%.2f',self.PSI),...
                        '\n\n']);
                end
                
                % generate (static) partership network structure
                    if self.VERBOSE
                        if isempty(adj_set)
                            fprintf(['Generating scale-free network structure with ALPHA = ',num2str(self.ALPHA),' ....']);
                        else
                            fprintf(['Using pre-generated network with ALPHA = ', num2str(adj_set.alpha_out),'....']);
                        end
                    end

                    tic;
                    if isempty(adj_set)
                        self.adj_full = [];
                        % loop mitigates minor bug producing errors for very small N<10 networks
                        while length(self.adj_full) ~= self.N   
                            [self.adj_full,~,self.rel_list_full] = self.pl_network(self.N,self.ALPHA,self.FULL_MAX_PARTNERS);
                        end

                    else
                        self.adj_full = adj_set.adj;
                        self.rel_list_full = adj_set.rel_list;
                        self.ALPHA = adj_set.alpha_out;
                        self.FULL_MAX_PARTNERS = adj_set.d_max;
                        self.adj_xmin = adj_set.x_min;
                        self.adj_L = adj_set.L;
                    end
                    network_elapsed = toc;
                

                % degree sequence of full network generated 
                % (from row sums of adjacency matrix)
                self.num_partners_full = full(sum(self.adj_full,2));

                % mean node degree of network
                self.mean_degree_full = mean(self.num_partners_full);

                % number, size and members of connected components
                [self.n_Comp,self.comp_sizes,self.comp_members] = self.networkComponents(self.adj_full);
                
                % set default daily partnership network
                self.adj = self.adj_full;
                self.num_partners = self.num_partners_full;
                
                % restrict daily partnership network by limiting -----
                % individuals to a max number of partners per day
                partner_check = self.num_partners_full - self.RESTRICT_MAX_PARTNERS;
                self.id_casuals = [find(partner_check > 0), partner_check(partner_check>0)];
                
                % pre-calculate fixed and prunable links
                self.fixed_rels = self.rel_list_full;
                
                % fixed links  (won't include any links to nodes with degree > threshold)
                for i = self.id_casuals(:,1)'
                    id_remove = self.fixed_rels(:,1) == i | self.fixed_rels(:,2) == i;
                    self.fixed_rels(id_remove,:) = []; 
                end
                % prunable links
                self.other_rels = setdiff(self.rel_list_full, self.fixed_rels,'rows');
               
                self.adj_int = sparse(self.N,self.N);
                self.mean_degree_current = zeros(1,2);
                self.restrict_adj3();
                %------------------------------------------------------
                
                % console output
                if self.VERBOSE
                    fprintf(['[DONE - in ',sprintf('%.2f',network_elapsed),' sec]']);
                    fprintf(['\n\t> mean degree (ave. partners) \t= ',sprintf('%.2f',self.mean_degree_full)]);
                    fprintf(['\n\t> # connected components \t= ',num2str(self.n_Comp),'\n\n']);
                end
                
              %% Infect population on day zero
                
                % start with non-AMR.. infect p0(1) fraction of population
                id_toinfect = randperm(self.N,round(self.p0(1)*self.N)); %random set of ids?
                self.indiv_state(id_toinfect,1,1) = 1; %infect people with those ids with nonAMR
                
                % for everyone infected with nonAMR gonorrhea at this point, 
                % switch a fraction p0(2) of those to be AMR infected instead
                idx = randperm(size(id_toinfect,2),round(self.p0(2)*size(id_toinfect,2)));
                self.indiv_state(id_toinfect(idx),2) = 1;
                self.indiv_state(id_toinfect(idx),1) = 0;
         
                % here take a proportion p0(3) of the AMR (only) infections
                % and (re)infect with nonAMR giving the prescribed
                % coinfection (given AMR)  proportion
                if self.ALLOW_COINFECTION == false && self.p0(3)~=0
                    self.p0(3) = 0;
                    warning('ALLOW_COINFECTION = false, setting p0(3) = 0');
                end
                
                id_AMR = find(self.indiv_state(:,2));
                idx = randperm(size(id_AMR,1),round(self.p0(3)*size(id_AMR,1)));
                self.indiv_state(id_AMR(idx),1) = 1;
                
                % CARE! what is the proportion of individuals who are
                % coinfected?? - how much difference does it make?
                %----------
                
                % set values for 'infected_since' and 'infection_count'
                % for each strain
                self.infected_since(self.indiv_state(:,1,1)==1,1) = 0;
                self.infected_since(self.indiv_state(:,2,1)==1,2) = 0;
                self.counters.infection_count(self.indiv_state(:,1,1)==1,1) = 1;
                self.counters.infection_count(self.indiv_state(:,2,1)==1,2) = 1;
                
                % prescribe proportion of symptomatic infections
                %size(self.indiv_state,2)
                idx_infected = any(self.indiv_state,2); %dim2 is strain so if there is a 1 in either strain
                size(idx_infected)
                self.symptoms(idx_infected,1) = rand(sum(idx_infected),1)<self.P_SYMPTOMS;
                
                % for initial infected, symptomatic individuals who seek
                % treatment - randomise treatment seeking day to be within
                % first month (30 days)
                idx_willseek = false(self.N,1);
                idx_willseek(logical(self.symptoms)) = rand(sum(self.symptoms),1)<self.P_SEEKS_TREATMENT;
                %self.seeks_treatment_on(idx_willseek) = randi([1,21],sum(idx_willseek),1);
                                
                % init. other counters for day zero
                self.counters.prevalence(1,:) = sum(self.indiv_state(:,:,1),1);
                self.counters.prev_both(1) = sum(all(self.indiv_state(:,:,1),2),1);
                self.counters.prev_either(1) = sum(any(self.indiv_state(:,:,1),2),1);
                self.counters.prev_nonAMRonly(1) = sum(self.indiv_state(:,1,1) & ~self.indiv_state(:,2,1),1);
                %self.counters.n_screened(1,:) = 0;
                %self.counters.n_traced(1,:) = 0;
                
                 % infected partner ratio (per degree)
                [self.ipr,~] = self.infected_partner_ratio(self.indiv_state(:,:,1),self.adj);

                
                % Population data updated for day zero
                self.today = 0;       
                
                self.burn_in = params.burn_in;
                
            end % class constructor
            
            
            function [self] = simulate(self,n_Days)
              %% Continue SIS simulation from where we left off
                % each iteration of the main loop computes the infection
                % state, and other counters, for the next day
                
                % resize arrays for new data keeping historical data up total #days if required
                if ~self.LOW_MEM
                    self.indiv_state = cat(3,self.indiv_state,zeros(self.N, self.n_Strains ,n_Days));
                end
                self.counters.prevalence = cat(1,self.counters.prevalence,zeros(n_Days,self.n_Strains));
                self.counters.prev_both = cat(1,self.counters.prev_both,zeros(n_Days,1));
                self.counters.prev_either = cat(1,self.counters.prev_either,zeros(n_Days,1));
                self.counters.prev_nonAMRonly = cat(1,self.counters.prev_nonAMRonly,zeros(n_Days,1));
                self.counters.incidence = cat(1,self.counters.incidence,zeros(n_Days,self.n_Strains));
                self.counters.incidence_either = cat(1,self.counters.incidence_either,zeros(n_Days,1));
                self.counters.incidence_new = cat(1,self.counters.incidence_new,zeros(n_Days,1)); 
                self.counters.cefta = cat(1,self.counters.cefta,zeros(n_Days,1));
                self.counters.cipr = cat(1,self.counters.cipr,zeros(n_Days,1));
                self.counters.cefta_notinf = cat(1,self.counters.cefta_notinf,zeros(n_Days,1));
                self.counters.cipr_notinf = cat(1,self.counters.cipr_notinf,zeros(n_Days,1));
                self.counters.cefta_AMR = cat(1,self.counters.cefta_AMR,zeros(n_Days,1));
                self.counters.cipr_AMR = cat(1,self.counters.cipr_AMR,zeros(n_Days,1));
                self.counters.cefta_nonAMR = cat(1,self.counters.cefta_nonAMR,zeros(n_Days,1));
                self.counters.cipr_nonAMR = cat(1,self.counters.cipr_nonAMR,zeros(n_Days,1));
                
                self.counters.n_attended_seeked = cat(1,self.counters.n_attended_seeked,zeros(n_Days,self.n_Strains+1));
                self.counters.n_attended_recalled = cat(1,self.counters.n_attended_recalled,zeros(n_Days,self.n_Strains+1));
                self.counters.n_attended_traced = cat(1,self.counters.n_attended_traced,zeros(n_Days,self.n_Strains+1));
                self.counters.n_attended_treated = cat(1,self.counters.n_attended_treated,zeros(n_Days,self.n_Strains+1));
                self.counters.n_treated_infected_symptomatic = cat(1,self.counters.n_treated_infected_symptomatic,zeros(n_Days,2));

                self.counters.n_screened = cat(1,self.counters.n_screened,zeros(n_Days,self.n_Strains+1));
                self.counters.n_traced = cat(1,self.counters.n_traced,zeros(n_Days,self.n_Strains+1));
                self.counters.n_recalled = cat(1,self.counters.n_recalled,zeros(n_Days,self.n_Strains+1));
                
                self.mean_degree_current = cat(1,self.mean_degree_current,zeros(n_Days,2));
                
                self.ipr = cat(1,self.ipr, nan(n_Days, max(self.num_partners)));

                %diary off;
                
                % counter only for this run of the simulation loop
                day_count = 0;
                
                % for console output only
                reverseStr = '';
                
                % start timer
                tic;
                
                % set defaults
                start_day = self.today;
                
                % day zero (self.today = 0) is at array position 1
                %for day = self.today+1:self.today+n_Days
                while self.today <= start_day + n_Days - 1
                    
                    % day counter for this simulation
                    day_count = day_count + 1;
                    self.today = self.today + 1;
                    
                    if self.DEBUG
                        fprintf(['\n-----------------------------\n','Day: ', num2str(self.today),'\n']);
                    end
                    
                    % update parameters if required (for changing model)
                    % dynamics after specified number of days
                    self.update_params();
                    
                    % update today's partnership network ever RESTRICT_RATE days
                    if mod(self.today-1,self.RESTRICT_RATE) == 0
                        self.restrict_adj3();
                    end
                    
                    % store current infection state and use as template for update state 
                    if self.LOW_MEM
                        current_state = self.indiv_state;
                    else                   
                        current_state = self.indiv_state(:,:,self.today);
                    end
                    
                    % set current state template which will be
                    % incrementally modified during each procedure below
                    new_state = current_state;


                  %% infection of strain susceptible individuals (zeros in state matrix)          
                    new_state = self.spread_infection(current_state, new_state);
                    
                  %%  Screen population - will not include infections acquired today
                    % update recall notifications for a proportion infected patients
                    
                    % screen entire population with rate / probability GAMMA
                    self.screening(current_state, true(self.N,1), self.GAMMA);  
                    

                  %% Parter-independent natural recovery with rate (r)
                    new_state = self.natural_recovery(current_state, new_state);
                    
                  %% Births / Deaths           
                    new_state = self.birth_death(new_state);

                  %% Treatment seeking
                    % 1) For symptomatic individuals seeking treatment, including 
                    %  previously symptomatic who have since recovered naturally
                    % : existing recalls / traces will also be checked
                    % 2) For recalled / traced individuals who have an
                    % 'appointment' for treatment today
                    if self.ALLOW_TREAT
                        new_state = self.seek_treatment(current_state, new_state);
                    end
                    
                  %% Finally: update state matrix for next day with 
                   % recoveries and infections produced by above procedures
                    if self.LOW_MEM
                        self.indiv_state = new_state;
                        % current prevalence (number of individuals infected with each strain)
                        % prevalences (as row vector)
                        self.counters.prevalence(self.today+1,:) = sum(self.indiv_state,1); %?? first dim is N?
                        self.counters.prev_both(self.today+1) = sum(all(self.indiv_state,2),1);
                        self.counters.prev_either(self.today+1) = sum(any(self.indiv_state,2),1);
                        self.counters.prev_nonAMRonly(self.today+1) = sum(self.indiv_state(:,1) & ~self.indiv_state(:,2),1);
                        
                        % update infected partner ratio
                        %[self.ipr(day+1,:),~] = self.infected_partner_ratio(self.indiv_state,self.adj);
                        
                    else
                        self.indiv_state(:,:,self.today+1) = new_state;
                        % current prevalence (number of individuals infected with each strain)
                        % prevalences (as row vector)
                        self.counters.prevalence(self.today+1,:) = sum(self.indiv_state(:,:,self.today+1),1);
                        self.counters.prev_both(self.today+1) = sum(all(self.indiv_state(:,:,self.today+1),2),1);
                        self.counters.prev_either(self.today+1) = sum(any(self.indiv_state(:,:,self.today+1),2),1);
                        self.counters.prev_nonAMRonly(self.today+1) = sum(self.indiv_state(:,1,self.today+1)...
                                                                      & ~self.indiv_state(:,2,self.today+1),1);
                        % update infected partner ratio
                        %[self.ipr(day+1,:),~] = self.infected_partner_ratio(self.indiv_state(:,:,day+1),self.adj);
                    end
                    
                    % clear symptoms of individuals who are no longer
                    % infected by any strain
                    idx_notinfected = ~any(new_state,2);
                    self.symptoms(idx_notinfected) = 0;
                    % other notifications (recalls/trace/seek) should
                    % have already been reset (if required) in the above procedures
                                       
                   

                  %% console output - show progress on screen
                    if self.VERBOSE && ~self.DEBUG
                        if self.burn_in
                            percentDone = day_count;
                            msg = [sprintf('Burning in... %3.0f', percentDone), ' days'];
                            fprintf([reverseStr, msg]);
                            reverseStr = repmat(sprintf('\b'), 1, length(msg));
                        else
                            percentDone = 100 * day_count / n_Days;
                            msg = [sprintf('Simulating... %3.1f', percentDone), '%%']; 
                            fprintf([reverseStr, msg]);
                            reverseStr = repmat(sprintf('\b'), 1, length(msg)-1);
                        end
                    end
                    
                  %% Burn-in (if selected)
                    % TESTING ONLY - NOT A STABLE METHOD OF EQUILIBRATING!!
                    if self.burn_in == true
                        % burn in
                        AMR_ratio = self.counters.prevalence(self.today,2) / self.counters.prev_either(self.today);

                        conditions = all([AMR_ratio >= 0.32]);
                        
                        if conditions == true || self.today == (start_day + n_Days - 1)
                           
                           % update state / counters for day zero with
                           % today's values
                           if self.LOW_MEM
                                self.indiv_state = new_state;
                                self.counters.prevalence(1,:) = sum(self.indiv_state,1);
                                self.counters.prev_both(1) = sum(all(self.indiv_state,2));
                                self.counters.prev_either(1) = sum(any(self.indiv_state,2));
                           else
                                self.indiv_state(:,:,1) = new_state;
                                self.counters.prevalence(1,:) = sum(self.indiv_state(:,:,self.today+1),1);
                                self.counters.prev_both(1) = sum(all(self.indiv_state(:,:,self.today+1),2));
                                self.counters.prev_either(1) = sum(any(self.indiv_state(:,:,self.today+1),2));
                           end
                           
                           % reset / update some totals
                           self.counters.drug_count = zeros(self.N,2);
                           self.infection_count = zeros(self.N,2);
                           self.counters.infection_count(self.indiv_state(:,1,1)==1,1) = 1;
                           self.counters.infection_count(self.indiv_state(:,2,1)==1,2) = 1;
                           
                           % adjust notifications / dates
                           self.infected_since = self.infected_since - self.today;
                           self.recall_notify(:,1:2) = self.recall_notify(:,1:2) - self.today;
                           self.trace_notify(:,1:2) = self.trace_notify(:,1:2) - self.today;
                           self.seeks_treatment_on = self.seeks_treatment_on - self.today;

                           
                           % reset day counters for this loop
                           start_day = 0;
                           self.today = 0;
                           
                           % burn-in completed?
                           if conditions == true
                               day_count = 0;
                               self.burn_in = false;
                           end
                           
                        end
                    end
                        

                end    
                loop_elapsed = toc;
                %diary on;

                % Simulation loop complete - output to console
                if self.VERBOSE; fprintf([' [DONE - in ',sprintf('%.2f',loop_elapsed),' sec]\n']);end
                 
                
            end % simulate function
            
            function [self] = update_params(self)
              %% Update current simulation parameters if required for studying model
                % dynamics under changing conditions
                % (function is called at the start of each new day)
                
                %ACTIVATION_DAY = 2*365;   % day on which parameters will change

                n_changes = size(self.param_updates,1);
                
                for k = 1:n_changes
                    
                    if self.today == self.param_updates{k,2}
                        evalc(self.param_updates{k,1});
                    end
                                     
                end              
                
            end
            
            function [new_state] = spread_infection(self, current_state, new_state)
              %% Infect population with each strain according to transmission rate 
                % and sate of the partnership network
                
                % NEED TO ADD:
                % day when infection acquired and whether it is symptomatic
                % with care not to overwrite existing data for infected
                % individuals!!!! 
%                 
                % DEBUG ONLY
%                 self = ibm;
%                 current_state = self.indiv_state(:,:,end);
                %--------------------------
                
                % calc. number of neighbours of each indiv. having each strain
                % [WARNING - Will be heavy-duty matrix multiplication for large N!]
                if self.LOW_MEM
                    nb_totals = self.adj*double(self.indiv_state);
                else                   
                    nb_totals = self.adj*double(self.indiv_state(:,:,self.today)); 
                end
                
                % rate of infection proportional to number of infected neighbours (per strain)
                % masking with ~current_state to ensure infection passes
                % only to those susceptible to the specific strain
                % OR
                % if coinfection not allowed, only individuals currently
                % not infected by either strain are infected (equal chance
                % for either nonAMR/ AMR strain
                % CARE: is EQUAL chance a reasonable assumption or should
                % this probability be weighted according to the relative
                % infection force (if partners have both strains to pass on)?
                if self.ALLOW_COINFECTION
                    mask = ~current_state;
                else
                    idx_not_infected = ~any(current_state,2);
                    mask = zeros(self.N,2);
                    mask(idx_not_infected,1) = rand(sum(idx_not_infected),1) < 0.5;
                    mask(idx_not_infected,2) = ~mask(idx_not_infected,1);
                end
                
                % Force of infection given number of partners
                % with each strain type
                infect_force = (1 - ( ( 1-repmat(self.BETA,self.N,1) ).^nb_totals )) .* mask;
                

                % susceptible individuals to infect with specific strain today
                idx_infect = rand(self.N,self.n_Strains) < infect_force;
                %idx_infect = sprand(infect_force) > (1-infect_force);


                % update temporary state matrix with infections
                % will be updated after other recovery processes completed
                % with history stored depending on global option LOW_MEM
                new_state(idx_infect) = 1;
                
                % update day of this infection (no history kept for now)
                self.infected_since(idx_infect) = self.today;
                
                % update infection count for each individual
                self.counters.infection_count(idx_infect) = self.counters.infection_count(idx_infect) + 1;
                
                % indices of individuals infected with either strain today
                idx_infected_today = any(idx_infect,2);
                                
                % assign symptoms according to probability of infected individual showing symptoms
                self.symptoms(idx_infected_today,1) = rand(sum(idx_infected_today,1),1) < self.P_SYMPTOMS;
                
                % indices of symptomatic individuals, infected today
                % (either strain) who will seek treatment (excluding those
                % already due to seek treatment, i.e. from infection with different strain)
                idx_willseek = false(self.N,1);
                %idx_willseek = idx_infected_today & self.symptoms...
                %                & isnan(self.seeks_treatment_on) & (rand(self.N,1) < self.P_SEEKS_TREATMENT);
                idx_willseek(idx_infected_today) = self.symptoms(idx_infected_today)...
                                & isnan(self.seeks_treatment_on(idx_infected_today))...
                                & (rand(sum(idx_infected_today),1) < self.P_SEEKS_TREATMENT);
                n_willseek = sum(idx_willseek,1);
                
                % normally distributed individual lab test and recall delays (integer rounded)
                indiv_onset_delays = round(normrnd(self.ONSET_DELAY_MEAN, self.ONSET_DELAY_STD,[n_willseek 1]));
                indiv_seek_delays = round(normrnd(self.SEEK_DELAY_MEAN, self.SEEK_DELAY_STD,[n_willseek 1]));

                % truncate delays at zero
                indiv_onset_delays(indiv_onset_delays < 0) = 0;
                indiv_seek_delays(indiv_seek_delays < 0) = 0;
                
                % add date of (voluntary / symptomatic) treatment seeking
                % to notification array
                self.seeks_treatment_on(idx_willseek) = self.today + indiv_onset_delays + indiv_seek_delays; 
                
                % update incidence of each strain for today 
                % (remember day zero is in array position 1
                    % new infections by each strain (columns)    
                     self.counters.incidence(self.today+1,:) = sum(idx_infect,1); 
                    % new infection by either strain
                     self.counters.incidence_either(self.today+1,:) = sum(any(idx_infect,2),1); 
                    % new infection by either strain excluding those already infected by one strain
                    % (i.e. new infection for an invi who was prev. susceptible to all strains)
                     self.counters.incidence_new(self.today+1,:) = sum(~any(current_state,2) & any(new_state,2),1); 
                
%                 if self.VERBOSE
%                     fprintf(['\nDay ',num2str(day),':\n',....
%                         '\t\t new non-AMR infections (', num2str(sum(idx_infect(:,1))),') :', num2str(find(idx_infect(:,1))),...
%                         '\n\t\t new AMR infections (', num2str(sum(idx_infect(:,2))),') :', num2str(find(idx_infect(:,2))),'\n\n']);
%                 end
                if self.DEBUG
                    fprintf(['\tNew infections: ', num2str(sum(idx_infect(:,1))),' (nonAMR), ',...
                        num2str(sum(idx_infect(:,2))),' (AMR)\n'])
                end
                
            end
            
            function [new_state] = natural_recovery(self, current_state, new_state)
              %% natural recovery of infected individuals AMR/non-AMR/coinfected
                % CARE: natural recovery is currently implemented for
                % strains independently (i.e. each strain has equal and 
                % independent prob. of recovering on any day

                % (excludes indivs. infected this time step)
                recover_rate = self.R*current_state;

                % individuals to recover today
                idx_recover = rand(self.N,self.n_Strains) < recover_rate;

                % update temproary state matrix with recoveries
                new_state(idx_recover) = 0;
                
                % individual now susceptible so clear 'infected since day' counted
                % (independently for each strain)
                self.infected_since(idx_recover) = NaN;
                
                % NOTE: notification flags (trace/recall/voluntary
                % treatment seeking) are all left unchanged
                
                if self.DEBUG
                    fprintf(['\tNatural recoveries: ', num2str(sum(idx_recover(:,1))),...
                        ' (nonAMR), ',num2str(sum(idx_recover(:,2))),' (AMR)\n']);
                end
                
            
            end
            
            function [new_state] = birth_death(self, new_state)
              %% recycle population with births and deaths 
                % - returning small fraction of population from infected (any strain) to
                % susceptible on each day
            
                % unlike the natural recovery of a single strain, the
                % birth or death of an individual should reset the
                % entire 'compartment' to be susceptible to all strains

                % individuals to recover via birth/death turnover today
                idx_birthdeath = rand(self.N,1) < self.MU;

                % update temproary state matrix with births/deaths
                % i.e. population renewal, infected --> susceptible (all strains)
                new_state(idx_birthdeath,:) = 0;
                
                % clear notifications for deceased individuals
                self.recall_notify(idx_birthdeath,:) = NaN;
                self.trace_notify(idx_birthdeath,:) = NaN;
                self.seeks_treatment_on(idx_birthdeath) = NaN;
                
                % replaced individual now susceptible so clear 'infected since day'counted
                self.infected_since(idx_birthdeath,:) = NaN;
                
                % clear drug counters for replaced individual
                self.counters.drug_count(idx_birthdeath,:) = 0;
                
                if self.DEBUG
                    fprintf(['\tBirths/deaths: ', num2str(sum(idx_birthdeath)),...
                        ' (recovered)\n']);
                end
            
            end
             
            function [self] = screening(self, current_state, idx_toscreen, screen_rate)
           %%  Screen population for infection
            % a fixed proportion of all individuals is screened 
                % DEBUG ONLY
%                 self = ibm;
%                 current_state = self.indiv_state(:,1:2,1);
%                 self.GAMMA = 0.3;
%                 self.LAB_DELAY = 7;
%                 self.RECALL_DELAY_MEAN = 7;
%                 self.RECALL_DELAY_STD = 1;

                % ------------
                
                % here 'idx_toscreen' will either be the whole population
                % i.e. = true(N,1) when performing daily screening with a
                % given rate,
                % OR  equal to a subset of the population, for example when
                % screening individuals who have been traced.
                
                % individuals to be screened today
                if screen_rate < 1
                    % screen individuals with given probability 
                    % (rate per day)
                    idx_screened = idx_toscreen & (rand(self.N,1) < screen_rate);
                else
                    % screen all selected individuals
                    idx_screened = idx_toscreen;
                end
               
                % update screened today counter
                self.counters.n_screened(self.today+1,:) = self.counters.n_screened(self.today+1,:) +...
                                                    [sum(idx_screened), sum(current_state(idx_screened,:),1)];

                
                % if screened AND infected (either strain) 
                % - individual will be notified (recalled) 
                idx_recalled = idx_screened & any(current_state,2);
                       
                % total number of individuals being recalled today
                n_recalled = sum(idx_recalled,1);
                
                % update counter for numbers recalled today (actual recall
                % date for individual will vary according to delays)
                self.counters.n_recalled(self.today+1,:) = self.counters.n_recalled(self.today+1,:) +...
                                                    [n_recalled, sum(current_state(idx_recalled,1)),...
                                                    sum(current_state(idx_recalled,2))];
                
                if n_recalled > 0
                    % some individuals may have an existing recall (having
                    % been screened again prior to attending for treatment)
                    idx_existing_recall = idx_recalled & ~isnan(self.recall_notify(:,1));

                    % new recalls
                    idx_new_recall = idx_recalled & ~idx_existing_recall;
                    n_new_recall = sum(idx_new_recall);

                    % normally distributed individual lab test and recall delays (integer rounded)
                    indiv_lab_delays = round(normrnd(self.LAB_DELAY_MEAN, self.LAB_DELAY_STD,[n_new_recall 1]));
                    indiv_recall_delays = round(normrnd(self.RECALL_DELAY_MEAN, self.RECALL_DELAY_STD,[n_new_recall 1]));

                    % truncate lab delays at 0
                    indiv_lab_delays(indiv_lab_delays < 0) = 0;
                    
                    % and recall delays at 1
                    indiv_recall_delays(indiv_recall_delays < 1) = 1;

                    % update notification array for new recalls with:
                    % 1. the day notified (today plus time taken for lab results)
                    % 2. day on which patient will return for treatment (truncated normal dist. above)
                    % 3. a flag indicating AMR (used actual AMR state from today)
                    self.recall_notify(idx_new_recall,:) = [self.today+indiv_lab_delays,...
                                                        self.today+indiv_lab_delays+indiv_recall_delays,...
                                                        current_state(idx_new_recall,2)];

%                     % for individuals with previous recall notice, only update
%                     % the AMR flag with latest results 
%                     % (use AMR strain flags directly)
                     self.recall_notify(idx_existing_recall,3) = current_state(idx_existing_recall,2);
%                     % NOTE: the above may change previously AMR flagged
%                     % notifications to nonAMR,  acceptable since this is based on a
%                     % more recent screening result
                end     
                
                % DEBUG CHECKS - run from console
                % histogram should give distribution of individual delay
                % times as prescribed by RECALL_DELAY_MEAN and RECALL_DELAY_STD
                
                if self.DEBUG
                    fprintf(['\tScreened: ',num2str(sum(idx_screened)),' (of which ',num2str(sum(idx_recalled)),' infected)\n']);
                    if n_recalled > 0
                        fprintf(['\t\t',num2str(n_new_recall),' (new recalls)\n']);
                        fprintf(['\t\t',num2str(n_recalled - n_new_recall),' (already recalled)\n']);
                    end
                end
                %check = ibm.recall_notify(ibm.recall_notify(:,1)>0,2)-ibm.recall_notify(ibm.recall_notify(:,1)>0,1);
                %hist(check);
                %------------------------
                
                                                
                
            %%
            end
            
            function [new_state] = seek_treatment(self, current_state, new_state)
              %%  treat symptomatic individuals who voluntarily report
                    % If an individual has been (recently) infected,
                    % reported symptoms and choses to seek treatment - they
                    % will have a 'date' flagged for treatment seeking
                    % (self.seeks_treatment_on)
                    % Here we check to see if today is their treatment day,
                    % treat accordingly, perform partner tracing, and reset
                    % any notifications
     
                    % DEBUG ONLY
%                      self = ibm;
%                      self.today = 13;
%                      current_state = self.indiv_state;
%                      new_state = self.indiv_state;
                    %-----------
                    
                    % which individuals are:
                    % symptomatics / prev. symptomatics seeking treatment today
                    idx_seeked = (self.seeks_treatment_on == self.today);
                    % were recalled and 'arrive' today
                    idx_recalled = (self.recall_notify(:,2) == self.today);
                    % were traced and 'arrive' today
                    idx_traced = (self.trace_notify(:,2) == self.today);
                    
                    % Default option: any individual who has voluntarily seeked
                    % treatment, been recalled or traced will be
                    % treated today (infected or not)
                    % [unless precreening or point-of-care enabled below]
                    idx_treat_today = any([idx_seeked, idx_recalled, idx_traced], 2);
                    
                    
                    % Option to screen traced individuals rather than treat
                    % today (ignore if we are using a point-of-care test)
                    idx_prescreen = false(self.N,1);              
                    if (self.PRESCREEN_TRACED || self.PRESCREEN_SEEKED) && ~self.ENABLE_POCT
                        % Rather than treat individuals who have been
                        % traced / seeked today (and not recalled today),
                        % screen (test) them first after which a recall 
                        % will be issued if found to be infected
                        
                        if self.PRESCREEN_TRACED && ~self.PRESCREEN_SEEKED
                            idx_prescreen = idx_traced & ~(idx_seeked | idx_recalled);
                        
                            % those remaining will be treated today 
                            idx_treat_today = any([idx_seeked, idx_recalled], 2);
                            
                            % clear these trace notifications for those who arrived
                            % today (thos infected will end up with a recall
                            % notification post-screening)
                            self.trace_notify(idx_prescreen,:) = NaN;
                            
                        elseif self.PRESCREEN_SEEKED && ~self.PRESCREEN_TRACED
                            idx_prescreen = idx_seeked & ~(idx_traced | idx_recalled);
                        
                            % those remaining will be treated today 
                            idx_treat_today = any([idx_traced, idx_recalled], 2);
                            
                            % clear these seek notifications for those who arrived
                            % today (thos infected will end up with a recall
                            % notification post-screening)
                            self.seeks_treatment_on(idx_prescreen) = NaN;
                            
                        elseif self.PRESCREEN_TRACED && self.PRESCREEN_SEEKED
                            idx_prescreen = idx_traced | idx_seeked & ~(idx_recalled);
                        
                            % those remaining will be treated today 
                            idx_treat_today = idx_recalled;
                            
                            % clear these trace notifications for those who arrived
                            % today (thos infected will end up with a recall
                            % notification post-screening)
                            self.trace_notify(idx_prescreen,:) = NaN;
                            self.seeks_treatment_on(idx_prescreen) = NaN;

                        end
                            
                        
                        % perform screening (assign recalls where infected)
                        self.screening(current_state, idx_prescreen, 1);                

                        % update screened today counter
                        self.counters.n_screened(self.today+1,:) = self.counters.n_screened(self.today+1,:) +...
                                                    [sum(idx_prescreen), sum(current_state(idx_prescreen,:),1)];
                        
                                              
                    end
                    
                    % number of individuals who will be treated today                    
                    n_treat_today = sum(idx_treat_today,1);
                    
                    % individuals who are treated and who are infected
                    % (either strain) - ONLY these indivs will be asked to
                    % provide partner details for tracing
                    idx_treated_infected = idx_treat_today & any(current_state,2);
                    n_treated_infected = sum(idx_treated_infected,1);
                    
                    % update counters for treated indivs. who were infected
                    % and how many were symptomatic
                    idx_treated_infected_symptomatic = idx_treated_infected & self.symptoms;
                    self.counters.n_treated_infected_symptomatic(self.today+1,:) = [sum(idx_treated_infected),...
                                                                        sum(idx_treated_infected_symptomatic)];
                    
                    % update some counters
                    self.counters.n_attended_seeked(self.today+1,:) = [sum(idx_seeked),...
                                                            sum(current_state(idx_seeked,:),1)];
                    self.counters.n_attended_recalled(self.today+1,:) = [sum(idx_recalled),...
                                                            sum(current_state(idx_recalled,:),1)];
                    self.counters.n_attended_traced(self.today+1,:) = [sum(idx_traced),...
                                                            sum(current_state(idx_traced,:),1)];          
                    %fprintf(num2str(self.counters.n_attended_traced(self.today+1,1)));
                    self.counters.n_attended_treated(self.today+1,:) = [sum(idx_treat_today),...
                                                            sum(current_state(idx_treat_today,:),1)];

                    % Proceed only if there are individuals seeking
                    % treatment today
                    if n_treat_today > 0          
                        
                        % compute who will be treated as non-AMR / AMR
                        % : start with a proprotion of cases in which AMR
                        % is assumed and (Ceft/A) is prescribed without 
                        % prior knowledge (currently around 86.5% of cases
                        % - given by self.P_BLINDTREAT_AS_AMR)
                        % This treatment yields (any) state [? ?] - [0 0]
                        % with the remainder treated with non-AMR (Cipr)
                        % yielding [? ?] --> [0 ?]
                        idx_treat_as_AMR = false(self.N,1);
                        idx_treat_as_nonAMR = false(self.N,1);
                        idx_treat_as_AMR(idx_treat_today) = rand(n_treat_today,1) < self.P_BLINDTREAT_AS_AMR;
                        % if proportion of cases with unknown AMR status
                        % treated as if AMR is very high (~1) then the
                        % list below will be empty
                        idx_treat_as_nonAMR(idx_treat_today) = ~idx_treat_as_AMR(idx_treat_today);
                       
                        allow_nonAMR = []; % select none by default
                        if (self.ENABLE_nonAMR_RECALL || self.ENABLE_nonAMR_TRACE) && ~self.ENABLE_POCT
                        % Determine if there are Recall / Trace notifications zeros
                        % which are relavent TODAY
                        %------------------------------------
                        % check if these individuals seeking treatement
                        % due to symptoms also have existing recall/trace
                        % notifications (treatment seeking individuals reporting
                        % due to symptoms may have already been recalled or
                        % traced).
                        % has_recall and has_trace below therefore list all
                        % individuals where recalls/traces apply
                            has_recall = false(self.N,1);
                            if self.ENABLE_nonAMR_RECALL
                                has_recall(idx_treat_today) = idx_recalled(idx_treat_today)...
                                    | ( idx_seeked(idx_treat_today)...
                                    & (self.recall_notify(idx_treat_today,1) < self.today) );
                            end

                            has_trace = false(self.N,1);
                            if self.ENABLE_nonAMR_TRACE
                                has_trace(idx_treat_today) = idx_traced(idx_treat_today)...
                                    | ( idx_seeked(idx_treat_today)...
                                    & (self.trace_notify(idx_treat_today,1) < self.today) );
                            end
                                % note that where a notification does not exist,
                                % values will be NaN and therefore fail the test (< self.today)

                            % If AMR notification status is flagged
                            % in EITHER notificiation then ensure always 
                            % treated with Ceft/A (and not non-AMR drug Cipr)
                            has_AMR_either = (has_recall & (self.recall_notify(:,3) == 1)... 
                                    | (has_trace & self.trace_notify(:,3) == 1));
                            idx_treat_as_AMR(has_AMR_either) = 1;
                            idx_treat_as_nonAMR(has_AMR_either) = 0;

                            % if however they were recalled or traced with non-AMR
                            % flag then we can make educated guess following
                            % some conditions, to prevent CEFT/A wastage.
                            has_nonAMR_both = (has_recall | has_trace) & ~has_AMR_either;
                            allow_nonAMR = find(has_nonAMR_both); 

                            % - CONDITIONS FOR extra nonAMR treatment pathway
                            % -------------------
                            if ~isempty(allow_nonAMR) 
                                % if difference between most recent recall/trace issue day and today
                                % exceeds threshold then will treat as AMR
                                delay_check = min( (self.today-[self.recall_notify(allow_nonAMR,1),...
                                                self.trace_notify(allow_nonAMR,1)]),[],2 ) > self.NON_AMR_MAX_DELAY;
                                % if individual has more than a specified
                                % number of partners then treat as AMR
                                partner_check = self.num_partners(allow_nonAMR) > self.NON_AMR_MAX_PARTNERS;

                                % remove entries which apply to the above
                                % (treating as if AMR)
                                allow_nonAMR(delay_check | partner_check) = [];

                                % what remains (if anything) is allowed to be treated as nonAMR
                                idx_treat_as_nonAMR(allow_nonAMR) = 1;
                                idx_treat_as_AMR(allow_nonAMR) = 0;
                            end
                        end
                        
                        % clear seek/recall/trace notifications
                        % for all individuals treated today
                        % (there may not be any at this point)
                        self.seeks_treatment_on(idx_treat_today) = NaN;
                        self.recall_notify(idx_treat_today,:) = NaN;
                        self.trace_notify(idx_treat_today,:) = NaN;
                        
                        
                        % Point-of care testing is enabled (treatment seekers
                        % today will be treated according to their infection
                        % status)
                        if self.ENABLE_POCT
                            if self.DISCRIM_POCT % Discriminatory POCT
                                % any with AMR strain
                                idx_treat_as_AMR = idx_treat_today & current_state(:,2); 
                            
                                % any with nonAMR strain ONLY
                                idx_treat_as_nonAMR = idx_treat_today & (current_state(:,1) & ~current_state(:,2));
                                
                            else % non-discriminatory POCT
                                
                                % treat any infected with either strain as AMR
                                idx_treat_as_AMR = idx_treat_today & any(current_state,2);
                                
                                % no nonAMR treatments
                                idx_treat_as_nonAMR = false(self.N, 1);
                                
                            end
                        end

                        % individuals treated with AMR therapy
                        %--------------------------------------------------------
                        
                        % set new state values both to zero (nonAMR and AMR)
                        % (cured/susceptible)
                        new_state(idx_treat_as_AMR,:) = 0;
                        
                        % clear 'infected since' log (both strains)
                        self.infected_since(idx_treat_as_AMR,:) = NaN;
                        
                        % increment CEFT/A counter for each individual
                        % treated today
                        self.counters.drug_count(idx_treat_as_AMR,2) = self.counters.drug_count(idx_treat_as_AMR,2) + 1;
                        
                        % log total # CEFT/A prescriptions for today
                        self.counters.cefta(self.today+1) = sum(idx_treat_as_AMR,1);
                        
                        % log number of 'correct' CEFT/A prescriptions
                        % (i.e. individuals had AMR infection)
                        self.counters.cefta_AMR(self.today+1) = sum(current_state(idx_treat_as_AMR,2) == 1, 1);
                        
                        % check and tally 'wasted' treatments - individuals
                        % : susceptibles [0 0] 
                        % (having recovered prior to treatment today)
                        self.counters.cefta_notinf(self.today+1) = sum(~any(current_state(idx_treat_as_AMR,:),2),1);
                        % : were non-AMR infected [1 0]
                        self.counters.cefta_nonAMR(self.today+1) = sum(current_state(idx_treat_as_AMR,1) == 1 &...
                                                                    current_state(idx_treat_as_AMR,2) == 0, 1);
                        % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                        
                        % individuals treated with non-AMR therapy
                        %-------------------------------------------------------------
                        % set non-AMR strain flags to zero
                        % (may still be AMR infected)
                        new_state(idx_treat_as_nonAMR, 1) = 0;

                        % increment cumulative drug (treatment) counters for each individual
                        self.counters.drug_count(idx_treat_as_nonAMR,1) = self.counters.drug_count(idx_treat_as_nonAMR,1) + 1;
                        
                        % log # Cipr prescriptions for today 
                        self.counters.cipr(self.today+1) = sum(idx_treat_as_nonAMR,1);
                        
                        % log number of 'correct' Cipr prescriptions
                        % (i.e. individuals had ONLY nonAMR infection)
                        self.counters.cipr_nonAMR(self.today+1) = sum(current_state(idx_treat_as_nonAMR,1) == 1 &...
                                                                    current_state(idx_treat_as_nonAMR,2) == 0, 1);
                        
                        % check and tally 'wasted' treatment
                        % : were susceptible [0 0]
                        % (having recovered prior to treatment)
                        self.counters.cipr_notinf(self.today+1) = sum(~any(current_state(idx_treat_as_nonAMR,:),2),1);
                        % : were AMR infected [? 1]
                        % (MISTREATED and will require patient recall)
                        idx_still_AMR = idx_treat_as_nonAMR & current_state(:,2)==1;
                        idx_correct = current_state(allow_nonAMR,1) & ~current_state(allow_nonAMR,2);
                        n_still_AMR = sum(idx_still_AMR);
                        self.counters.cipr_AMR(self.today+1) = n_still_AMR;
                        
                        % get extra info on these undertreated cases
                        %if n_still_AMR > 0
                            self.counters.informed.undertreated_degree = cat(1,self.counters.informed.undertreated_degree,full(sum(self.adj(idx_still_AMR,:),2)));
                            self.counters.informed.undertreated_infected_duration = cat(1,self.counters.informed.undertreated_infected_duration, self.today - self.infected_since(idx_still_AMR,:));
                            self.counters.informed.correct_degree = cat(1,self.counters.informed.correct_degree,full(sum(self.adj(idx_correct,:),2)));
                            self.counters.informed.correct_infected_duration = cat(1,self.counters.informed.correct_infected_duration, self.today - self.infected_since(idx_correct,:));
                            
                            %end
                         
                        % clear 'infected since' log (only nonAMR strain)
                        self.infected_since(idx_treat_as_nonAMR,1) = NaN;
                        
                        % normally distributed individual lab test delays 
                        % (integer rounded, truncated with a min. of 0 days)
                        % (compute value for every individual - bit wasteful)
                        indiv_lab_delays = round(normrnd(self.LAB_DELAY_MEAN, self.LAB_DELAY_STD,[self.N 1]));
                        indiv_lab_delays(indiv_lab_delays < 0) = 0;
                        
                        % RECALL PATIENTS who were mis-treated (had AMR!)
                        if n_still_AMR > 0
                            % normally distributed individual recall delays (integer rounded)
                            indiv_recall_delays = round(normrnd(self.RECALL_DELAY_MEAN, self.RECALL_DELAY_STD,[n_still_AMR 1]));

                            % truncate recall delays at 1
                            indiv_recall_delays(indiv_recall_delays < 1) = 1;

                            % RECALL notification array for (mistreated) AMR patients:
                            % 1. the day notified (today plus time taken for lab results (truncatd normal dist.)
                            % 2. day on which patient will return for treatment (truncated normal dists. above)
                            % 3. a flag indicating AMR (always 1 here!)
                            self.recall_notify(idx_still_AMR,:) = [self.today+indiv_lab_delays(idx_still_AMR),...
                                                                self.today+indiv_lab_delays(idx_still_AMR)+indiv_recall_delays,...
                                                                ones(n_still_AMR,1)];
                                                            
                            % update counter for numbers recalled today (actual recall
                            % date for individual will vary according to delays)
                            self.counters.n_recalled(self.today+1,:) = self.counters.n_recalled(self.today+1,:) +...
                                                    [n_still_AMR, sum(current_state(idx_still_AMR,1)),...
                                                    sum(current_state(idx_still_AMR,2))];
                                                            
%                             % If these individuals had a symptomatic
%                             % infection then it may continue to be
%                             % symptomatic after (incorrect treatment)
%                             % : set seek_treatment_on flags
%                             
%                             % normally distributed individual seek delays (integer rounded)
%                             indiv_seek_delays = round(normrnd(self.SEEK_DELAY_MEAN, self.SEEK_DELAY_STD,[n_still_AMR 1]));
% 
%                             % truncate seek delays at zero
%                             indiv_seek_delays(indiv_recall_delays < 0) = 0;
%                             
%                             % set 'seeks treatment on' flags (day to treat)
%                             % : note there is no additional symptom onset
%                             % delay here. If patient returns before lab
%                             % results (and therefore recall issued) - may
%                             % be treated again as nonAMR?!
%                             self.seeks_treatment_on(idx_still_AMR) = self.today + indiv_seek_delays;
                            
                            
                            % DEBUG:
                            %display([self.recall_notify(idx_still_AMR,:),current_state(idx_still_AMR,:)])
                                                            
                        end
                       
                        
                        % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                        
                        % TRACE PARTNERS - depending on AMR status
                        % want to trace partners of individuals treated
                        % today - who were infected (either strain)                  
                        % Tracing is not 100% efficient (not all partners
                        % of all individuals can be contacted)
                        % : create new 'partial' adjacency with only a
                        % fraction (self.PSI) of each individuals partners  
                        r = sprand( self.adj(:,idx_treated_infected) );
                        partial_adj = r > (1-self.PSI);
 
                        % Apply threshold self.MAX_TRACE which defines
                        % maximum number of individuals who can be traced
                        % from each index patient (randomised pruning where
                        % more than threshold have been selected so far)
                        
                        % values > 0 here indicate columns which need pruning
                        nb_sum_thresh = sum(partial_adj,1)-self.MAX_TRACE;

                        % get array subindices of of all entries in partial matrix
                        [i,j] = ind2sub(size(partial_adj),find(partial_adj));

                        % only prune for individuals  where threshold exceeded
                        for k = find(nb_sum_thresh>0)

                            % row indices (tracable partners) of each index patient where more
                            % pruning is needed
                            id_all = i(j==k);

                            % randomly selected a number of these according to how many need to be
                            % pruned
                            id_rmv = randperm(numel(id_all),nb_sum_thresh(k));

                            % remove selected entries from partial adjacency matrix
                            partial_adj(id_all(id_rmv),k) = 0;

                        end
%                         % values > 0 here indicate columns which need pruning
%                         check = sum(partial_adj,1)-self.MAX_TRACE;
%                         if any(check>0)
%                             warning('More than specified number of partners are being traced!!');
%                         end
                        
                        
                        % Nx2 matrix representing how many partners of each
                        % individual have been detected with each strain
                        trace_count = partial_adj * double(current_state(idx_treated_infected,:));
                                
                        % individuals who will be notified today
                        % CARE - probably should exclude everyone treated
                        % today?? Does it depend on AMR status?
                        idx_trace = any(trace_count,2);
                        n_to_trace = sum(idx_trace,1);
                        
                        % update trace today counter
                        self.counters.n_traced(self.today+1,:) = [sum(idx_trace),...
                                                                sum(current_state(idx_trace,:),1)];
                        
                        
                        if n_to_trace > 0
                            
                            % some individuals may already have a trace
                            % notification (these notifications will only
                            % modified if there a new AMR risk)
                            idx_existing_trace = idx_trace & ~isnan(self.trace_notify(:,1));
                            n_existing_trace = sum(idx_existing_trace);
                            idx_new_trace = idx_trace & ~idx_existing_trace;
                            n_new_trace = n_to_trace - n_existing_trace;
                            
                            if (self.LAB_DELAY_MEAN == 0) && (self.LAB_DELAY_STD ==0)
                                % LAB delay is zero for all traces
                                min_lab_delays = zeros(n_new_trace,1);
                            else
                                % some individuals may be traced by multiple
                                % infectees today - from each of the index cases
                                % find the minimum lab delay
                                delay_mat = full( self.adj(idx_new_trace,idx_treated_infected)...
                                    .*repmat(indiv_lab_delays(idx_treated_infected)',n_new_trace,1) );
                                delay_mat(delay_mat==0) = inf;
                                % final lab delay value for each traced individual
                                min_lab_delays = min(delay_mat,[],2);
                            end
                            
                            % individual trace delays (normal. dist. truncated at 1)
                            indiv_trace_delays = round(normrnd(self.TRACE_DELAY_MEAN, self.TRACE_DELAY_STD,[n_new_trace 1]));
                            indiv_trace_delays(indiv_trace_delays < 1) = 1;

                            % update TRACE notification array for new traces with:
                            % (excluding those will exising trace notifications)
                            % 1. the day notified (today plus time taken for lab results (truncatd normal dist.)
                            % 2. day on which patient will return for treatment (truncated normal dists. above)
                            % 3. a flag indicating whether trace is nonAMR (0) or  AMR (1)
                            %    determined by positive values in the 2nd column of trace_count  
                            self.trace_notify(idx_new_trace,:) = [self.today+min_lab_delays,...
                                        self.today+min_lab_delays+indiv_trace_delays,...
                                        trace_count(idx_new_trace,2)>0];
                                    
                            % for individuals with existing trace
                            % notifications - update those where any of the
                            % traces imply an AMR infection by modifying
                            % the AMR flag only (existing dates remain -
                            % trace may already have been AMR)
                            current_flag = self.trace_notify(idx_existing_trace,3);
                            new_flag = trace_count(idx_existing_trace,2)>0;
                            self.trace_notify(idx_existing_trace,3) = any([current_flag new_flag],2);
                            
                        end
                        
                    end % individuals seeking treatment today
                    
                    % Output extra info to console for debugging
                    if self.DEBUG
                        fprintf(['\tTreated today: ',num2str(sum(idx_treat_today)),'\n']);
                            if n_treat_today > 0
                                fprintf(['\t\t',num2str(sum(idx_seeked)),' (seeked voluntarily)\n']);
                                fprintf(['\t\t',num2str(sum(idx_recalled)),' (recalled)\n']);
                                fprintf(['\t\t',num2str(sum(idx_traced)),' (traced)\n']);
                                if self.PRESCREEN_TRACED
                                    fprintf(['\t\t\t(of which ',num2str(sum(idx_prescreen)),' screened instead of treated)\n']);
                                end
                                fprintf('\t\t---------------------\n');
                                fprintf(['\t\t',num2str(sum(idx_treat_as_nonAMR)),' (treated as if nonAMR)\n']);
                                fprintf(['\t\t',num2str(sum(idx_treat_as_AMR)),' (treated as if AMR)\n']);
                                fprintf(['\t\t',num2str(n_still_AMR),' recalled (wrong treatment)\n']);
                                if n_to_trace > 0
                                    fprintf('\t\t---------------------\n');
                                    fprintf(['\t\t',num2str(n_to_trace),' individuals traced\n']);
                                    fprintf(['\t\t\t',num2str(sum(n_new_trace)),' (newly traced indivs.)\n']);
                                    fprintf(['\t\t\t',num2str(sum(n_existing_trace)),' (existing traces updated)\n']);
                                end
                            end
                    end
                   
                
            end
            
%             function [ self ] = restrict_adj(self)
%           %% Returns a restricted partnership network where individiuals
%             % referenced in idx_mod have n_to_remove partnerships removed
%             % randomly from the original adjacency matrix
% 
%                 % check to see if we can skip this reduction 
%                 if self.RESTRICT_MAX_PARTNERS < self.FULL_MAX_PARTNERS
% 
% 
%                     % use full matrix representation for faster indexing
%                     % (will use more memory)
%                     adj_out = full(self.adj_full);
% 
%                     % get list from global value to start with
%                     self.id_casuals = self.id_casuals;
% 
%                     % LOOP START----
%                     while size(self.id_casuals,1) > 0
%                         % random entry from id_casuals
%                         r = randi(size(self.id_casuals,1),1);
% 
%                         this_node = self.id_casuals(r,1);
%                         this_nbs = find(adj_out(this_node,:));  % SLOWWWW
% 
%                         % randomly pick links to remove
%                         rm = randperm(numel(this_nbs), self.id_casuals(r,2));
%                         those_nodes = this_nbs(rm);
% 
%                         % add i,j th elements to set to zero this iteration 
%                         j = repmat(this_node,length(those_nodes),1);
%                         k = those_nodes';
% 
%                         % set elements (symmetrically) to zero   SLOWW!!!!
%                         adj_out(sub2ind([self.N,self.N],j,k)) = 0; 
%                         adj_out(sub2ind([self.N,self.N],k,j)) = 0;
% 
%                         % connections between 'this_node' and 'those_nodes' will be 
%                         % severed
%                         % i.e. degree of 'this_node' will reduce by number of 'those_nodes'
%                         % and therefore can be removed from list
%                         self.id_casuals(r,:) = [];
% 
%                         % degree of 'those_nodes' will subsequenly each be reduced by one
%                         % so adjust counters in id_casuals if necessary
%                         for l = 1:length(those_nodes)
%                             idx  = ( self.id_casuals(:,1)== those_nodes(l) );
%                             self.id_casuals(idx,2) = self.id_casuals(idx,2) - 1;
%                         end
% 
%                         % if these nodes now have degree less than the threshold,
%                         % we can also remove them from the list 'id_casuals'
%                         self.id_casuals(self.id_casuals(:,2)<=0,:) = [];
% 
%                     end
% 
%                     % degree distribution of reduced partnership network
%                     self.num_partners = sum(adj_out,2);
% 
%                     % assign 'todays' partnerhip network as the sparse ouput matrix
%                     self.adj = sparse(adj_out);
% 
% %                    %% checks
% %                     sym_out = issymmetric(adj_out)
% %                     max_deg_old = max(full(sum(self.adj_full,2)))
% %                     max_deg_new = max(full(sum(adj_out,2)))
% 
%                     if max(self.num_partners) > self.RESTRICT_MAX_PARTNERS
%                         error('Network restriction failed!');
%                     end
%                     
%                 end
%                 
%                 
% 
%                 
%             end
%             
%    
%             function [ self ] = restrict_adj2(self)
%           %% Returns a restricted partnership network where individiuals
%             % referenced in idx_mod have n_to_remove partnerships removed
%             % randomly from the original adjacency matrix
% 
%                 % check to see if we can skip this reduction 
%                 if self.RESTRICT_MAX_PARTNERS < self.FULL_MAX_PARTNERS
%                     
%                     % start will full undirected list of relationships
%                     new_rels = self.rel_list_full;
%     
%                     % loop through nodes which require link pruning
%                     for j = randperm(size(self.id_casuals,1))
%                         id = find(new_rels(:,1)==self.id_casuals(j,1) | new_rels(:,2)==self.id_casuals(j,1));
%                         d = size(id,1);
%                         rm = d-self.RESTRICT_MAX_PARTNERS;
%                         rm(rm<0) = 0;
%                         r = randperm(size(id,1),rm);
%                         new_rels(id(r),:) = [] ;   
%                     end
%                     
%                     % generate new undirected sparse matrix from pruned
%                     % relationship list
%                     self.adj = sparse( [new_rels(:,1);new_rels(:,2)], [new_rels(:,2);new_rels(:,1)], 1, self.N, self.N );
%                     
%                     % update the degree sequence (number of current
%                     % partners for each indiv)
%                     self.num_partners = full(sum(self.adj,2));
%                     
%                     % update mean degree history
%                     %self.mean_degree_current(self.today:self.today+self.RESTRICT_RATE-1) = mean(self.num_partners);
%                     self.mean_degree_current(self.today:end,:) = repmat([mean(self.num_partners),std(self.num_partners)],size(self.mean_degree_current(self.today:end,:),1),1);
% 
%                     % update the integrated network variable
%                     % (keeping track of all historical partners up until
%                     % now)
%                     self.adj_int = self.adj_int | self.adj;
%     
% 
% %                    %% checks
% %                     sym_out = issymmetric(adj_out)
% %                     max_deg_old = max(full(sum(self.adj_full,2)))
% %                     max_deg_new = max(full(sum(adj_out,2)))
%                     if max(self.num_partners) > self.RESTRICT_MAX_PARTNERS
%                         error('Network restriction failed!');
%                     end
%                 end
% 
%                 
%             end
            
            function [ self ] = restrict_adj3(self)
          %% Returns a restricted partnership network where individiuals
            % referenced in idx_mod have n_to_remove partnerships removed
            % randomly from the original adjacency matrix

                % check to see if we can skip this reduction 
                if self.RESTRICT_MAX_PARTNERS < self.FULL_MAX_PARTNERS
                    
                    % default - no links accepted
                    idx = false(length(self.other_rels),1);
                    
                    % zero degrees for (all) nodes
                    ds_check = zeros(self.N,1);
                    
                    % loop - check each link in random order
                    for i = randperm(length(self.other_rels))

                        % this relationship (two nodes)
                        rel = self.other_rels(i,:)';
                        
                        if max(ds_check(rel)) < self.RESTRICT_MAX_PARTNERS
                            % accept link if degrees are under threshold
                            idx(i) = 1;
                            
                            % increment degree of both nodes
                            ds_check(rel) = ds_check(rel)+1;
                        end  

                    end

                    % combine fixed edges with allowed random edges
                    new_rels = cat(1,self.fixed_rels,self.other_rels(idx,:));

                    % update adjacency matrix
                    self.adj = sparse( [new_rels(:,1); new_rels(:,2)], [new_rels(:,2); new_rels(:,1)], 1, self.N, self.N);

                    % update degree sequence of current network
                    self.num_partners = full(sum(self.adj,2));
                    
                    % update mean degree history
                    %self.mean_degree_current(self.today:self.today+self.RESTRICT_RATE-1) = mean(self.num_partners);
                    self.mean_degree_current(self.today:end,:) = repmat([mean(self.num_partners),std(self.num_partners)],size(self.mean_degree_current(self.today:end,:),1),1);
                    
                    % update the integrated network variable
                    % (keeping track of all historical partners up until
                    % now)
                    self.adj_int = self.adj_int | self.adj;

                    % checks
%                     sym_out = issymmetric(self.adj)
%                     max_deg_old = max(full(sum(self.adj_full,2)))
%                     max_deg_new = max(full(sum(self.adj,2)))

    
                    if max(self.num_partners) > self.RESTRICT_MAX_PARTNERS
                        error('Network restriction failed!');
                    end
                    
                end

            end
            
            function [self] = plots(self)
              %% Call static function plot_data to plot various outputs 
              % in 'counters' from day 0 until today
                
                if self.VERBOSE;fprintf('\nGenerating plots...');end
                
                % return if no data has been simulated
                if self.today < 1
                    if self.VERBOSE;fprintf('...[NO DATA]\n');end
                    return
                end
                
                ave = 7;
                
                [h_prev, n_prev] = self.plot_prev(self.counters, self.N, [0 self.today],self.DIR_NAME);
                [h_inc, n_inc] = self.plot_incidence(self.counters, self.N, [0 self.today], ave, self.DIR_NAME);
                [h_doses, n_doses] = self.plot_doses(self.counters, self.N, 7,self.DIR_NAME);
                [h_attend, n_attend] = self.plot_attendance(self.counters, self.N, [0 self.today], ave, self.DIR_NAME);
                [h_AMR_ratio, n_AMR_ratio] = self.plot_AMR_ratio(self.counters, self.DIR_NAME);
                
                self.fig_h = [h_prev, h_inc, h_doses, h_attend, h_AMR_ratio];
                self.plot_names = {n_prev, n_inc, n_doses, n_attend, n_AMR_ratio};
                
                if self.VERBOSE;fprintf('[DONE]\n');end           
              
            end
            
            function [self] = save_data(self)
            %% Save data files (all data) and plots if available
           
                % console output
                if self.VERBOSE;fprintf('\nSaving...');end
            
                % check output folders exist, if not then create
                if exist('./outdata','dir')~=7;mkdir('./outdata');end
                if exist(['./outdata/',self.DIR_NAME],'dir')~=7;mkdir(['./outdata/',self.DIR_NAME]);end
                
%                 % update file suffix
%                 self.DIR_NAME = ['_N=',num2str(self.N),'_nS=',num2str(self.n_Strains),...
%                                 '_days=',num2str(self.today)];
%                 
%                     for i = 1:length(self.fig_h)
%                         % if figure still alive then save
%                         if ishandle(self.fig_h(i))            
%                             % use export_fig for pdf output if available
%                             if exist('export_fig','file')~=0
%                                 export_fig(['./outdata/plots/',self.plot_name{i},self.DIR_NAME,'.pdf'],self.fig_h(i))
%                                 export_fig(['./outdata/plots/',self.plot_name{i},self.DIR_NAME,'.png'],'-r200',self.fig_h(i))
%                             end
%                             % also save matlab 'fig' file
%                             saveas(self.fig_h(i),['./outdata/plots/',self.plot_name{i},self.DIR_NAME,'.fig']);                            
%                         end  
%                     end
                

                % save simulation datafile (everything for now!);
                save(['./outdata/',self.DIR_NAME,'/sim_data.mat']);
                % and a copy of the class file
                copyfile('AMR_IBM.m',['./outdata/',self.DIR_NAME,'/AMR_IBM.m']);
                
                % console output
                if self.VERBOSE;fprintf('[DONE]\n\n');end
    
            end
            
            function [vals] = get_infected_per_degree(self)
              %% returns the number of individuals infected with each strain
                % or coinfecteds, as a function of node degree (number of partners)
                
                % OUTPUT:   vals    -  list of # infecteds with column for each 
                %                       strain with 3rd column for
                %                       coinfecteds. Row 1 = degree 1
                %                       indivs. etc.
                %                      (if there are no indivs. of degree
                %                      n, entire row will be NaNs)
                
                
                % DEBUG ONLY]
                self = ibm;
                
                %--------
              
                
                
                %max_d = self.FULL_MAX_PARTNERS;
                max_d = max(self.num_partners);
                vals = zeros(max_d,3);
                for d = 1:max_d
                    % infection state of all indivs. with degree 'd'
                    s = self.indiv_state(self.num_partners == d,:);
                    
                    if ~isempty(s)
                        % row for this degree containing 
                        % 1: # nonAMR, 2: # AMR, 3: # coinfectedthen
                        vals(d,:) = [sum(s,1) sum(all(s,2),1)];
                    else
                        % if no individuals of this degree then fill with
                        % NaNs
                        vals(d,:) = nan;
                    end
                end
                
            end
            
          
        end % class methods
        
        methods(Static)
            
            function [ adj, edge_list, rel_list ] = pl_network(N,alpha,d_max)
              %% Generate network with power-law degree distribution 
                %  such that the probability of a node having k neighbours 
                %  has the form P(k) - k^(-alpha)
                %
                %  Resulting network is undirected (symmetric) such that adj=adj'
                %  and has no self-loops (0 diagonal)
                % 
                %  INPUTS:
                %       N           - Number of nodes required.
                %       alpha       - power-law slope
                %       d_max       - maximum degree of any node (minimum = 1)
                %
                %  OUTPUTS:
                %       adj         - NxN network adjacency (connectivity) matrix
                %       edge_list   - full list of network edges (Nx3)
                %                       where 3rd column gives edge weight = 1
                %then
                %  BUGS:    has been found to occasionally produce adjacency matrix with one less
                %           row/colum than required by N for small N < 10. 
                %
                %  AUTHOR:  Adam Zienkiewicz (2016)  adamziek@gmail.com


              %% generate node degree distribution

                % min node degree
                d_min = 1;
 
                % power-law distributed sequence between upper and lower bounds
                % (see matlab)
                node_degree = ( ( d_max^(-alpha+1) - d_min^(-alpha+1) ).*rand(N,1) + d_min^(-alpha+1) ).^( 1/(-alpha+1) );

                % get integer node degree sequence (how accurate is rounding here?)
                % should truncation be avoided instead rounding to nearest
                % int?
                int_node_degree = fix(node_degree);
                
                % network constructor
   
                % init. list of free stubs for each node
                free_stubs = int_node_degree;

                % init edge list
                edge_list = [];
                rel_list = [];
                
                % sparse adjacency matrix
                adj = logical(sparse(N,N));

                % lookup vectors for each node's connections - too avoid self-loop and duplicates
                edge_lookup = num2cell((1:N)');

                % sequential loop through initial nodes
                % Does this need to be randomised? Since degree sequence is itself random, 
                % I *doubt* it makes any difference
                for this_node = randperm(N);%1:N   

                    % number of remaining stubs for selected node
                    my_stubs = free_stubs(this_node);

                    % free stubs available - excluding self and previously connected nodes
                    other_stubs = free_stubs;
                    other_stubs(edge_lookup{this_node}) = 0;

                    % loop over this node's free stubs
                    while my_stubs > 0
                        % indices of (other) nodes with free stubs
                        idx_avail_stubs = find(other_stubs > 0);

                        % if there is an available stub 
                        if ~isempty(idx_avail_stubs)
                            % chose one at random
                            that_node = idx_avail_stubs(randperm(length(idx_avail_stubs),1));
                            
                            % add to edge list (both directions)
                            edge_list = [edge_list;[this_node that_node];[that_node this_node]];
                            
                            % add to relationship list (one way)
                            rel_list = [rel_list;[this_node that_node]];
                            
                            % add to adjacency matrix
                            adj(this_node,that_node) = 1;
                            adj(that_node,this_node) = 1;

                            % reduce stub count for this node within this loop
                            my_stubs = my_stubs-1;

                            % reduce stub count for this node and connected node globally
                            free_stubs([this_node,that_node]) = free_stubs([this_node,that_node])-1;

                            % update lookup list - avoid duplicate edges
                            edge_lookup{this_node} = [edge_lookup{this_node} that_node];
                            edge_lookup{that_node} = [edge_lookup{that_node} this_node];

                            % update available stubs to this node
                            % free stubs available -  excluding this node and recently connected
                            other_stubs = free_stubs;
                            other_stubs(edge_lookup{this_node}) = 0;

                        else
                            % cannot find any free nodes for remaining stubs - break from this node
                            break
                        end

                    end % loop over node stubs


                end  % loop over selected node

                % add column of ones to edge list
                edge_list = [edge_list,ones(length(edge_list),1)];
                
                % tests
                [alpha_cont,x_cont,L_cont] = AMR_IBM.plfit(node_degree);
                [alpha_int,x_int,L_int] = AMR_IBM.plfit(int_node_degree);
                [alpha_adj,x_adj,L_adj] = AMR_IBM.plfit(full(sum(adj,2)));
                display(['alpha (continuous): ',num2str(alpha_cont)]);  
                display(['alpha (integers): ',num2str(alpha_int)]);
                display(['alpha (network): ',num2str(alpha_adj)]);
       
            end  % function
            
            function [ adj, edge_list, rel_list ] = pl_network2(N,alpha,d_max)
              %% Generate network with power-law degree distribution 
              
                %
                %  DO NOT USE THIS FUNCTION- TESTING ONLY!!!
                %
                
                %  such that the probability of a node having k neighbours 
                %  has the form P(k) - k^(-alpha)
                %
                %  Resulting network is undirected (symmetric) such that adj=adj'
                %  and has no self-loops (0 diagonal)
                % 
                %  INPUTS:
                %       N           - Number of nodes required.
                %       alpha       - power-law slope
                %       d_max       - maximum degree of any node (minimum = 1)
                %
                %  OUTPUTS:
                %       adj         - NxN network adjacency (connectivity) matrix
                %       edge_list   - full list of network edges (Nx3)
                %                       where 3rd column gives edge weight = 1
                %then
                %  BUGS:    has been found to occasionally produce adjacency matrix with one less
                %           row/colum than required by N for small N < 10. 
                %
                %  AUTHOR:  Adam Zienkiewicz (2016)  adamziek@gmail.com


              %% generate node degree distribution

                % min node degree
                d_min = 1;
 
                % power-law distributed sequence between upper and lower bounds
                % (see matlab)
                node_degree = ( ( d_max^(-alpha+1) - d_min^(-alpha+1) ).*rand(N,1) + d_min^(-alpha+1) ).^( 1/(-alpha+1) );

                % get integer node degree sequence (how accurate is rounding here?)
                % should truncation be avoided instead rounding to nearest
                % int?
                int_node_degree = fix(node_degree);
                
                % network constructor
   
                % init. list of free stubs for each node
                free_stubs = int_node_degree;

                % init edge list
                edge_list = [];
                rel_list = [];
                
                % sparse adjacency matrix
                adj = logical(sparse(N,N));

                % lookup vectors for each node's connections - too avoid self-loop and duplicates
                edge_lookup = num2cell((1:N)');
                
                % boolean list of complete nodes
                nodes_complete = free_stubs<=0;

                % cycle  through nodes until all are complete
                while ~all(nodes_complete)
                    node_list = find(nodes_complete==0);
                    this_node = node_list(randi(length(node_list),1));

                    % number of remaining stubs for selected node
                    my_stubs = free_stubs(this_node);

                    % free stubs available - excluding self and previously connected nodes
                    other_stubs = free_stubs;
                    other_stubs(edge_lookup{this_node}) = 0;

                    % loop over this node's free stubs
                    %while my_stubs > 0
                    % indices of (other) nodes with free stubs
                    idx_avail_stubs = find(other_stubs > 0);

                    % if there is an available stub 
                    if ~isempty(idx_avail_stubs)
                        % chose one at random
                        that_node = idx_avail_stubs(randperm(length(idx_avail_stubs),1));

                        % add to edge list (both directions)
                        edge_list = [edge_list;[this_node that_node];[that_node this_node]];

                        % add to relationship list (one way)
                        rel_list = [rel_list;[this_node that_node]];

                        % add to adjacency matrix
                        adj(this_node,that_node) = 1;
                        adj(that_node,this_node) = 1;

                        % reduce stub count for this node within this loop
                        my_stubs = my_stubs-1;

                        % reduce stub count for this node and connected node globally
                        free_stubs([this_node,that_node]) = free_stubs([this_node,that_node])-1;

                        % update lookup list - avoid duplicate edges
                        edge_lookup{this_node} = [edge_lookup{this_node} that_node];
                        edge_lookup{that_node} = [edge_lookup{that_node} this_node];

                        % update available stubs to this node
                        % free stubs available -  excluding this node and recently connected
                        other_stubs = free_stubs;
                        other_stubs(edge_lookup{this_node}) = 0;

                    else
                        % cannot find any free nodes for remaining stubs - break from this node
                        nodes_complete(this_node) = 1;
                        break
                    end
                end
                
                % add column of ones to edge list
                edge_list = [edge_list,ones(length(edge_list),1)];

       
            end  % function
            
            function [nComponents,sizes,members] = networkComponents(A)
              %% get number, size and membership of network connected components
                % Daniel Larremore
                % April 24, 2014
                % larremor@hsph.harvard.edu
                % http://danlarremore.com
                % Comments and suggestions always welcome.
                % 
                % INPUTS:
                % A                     Matrix. This function takes as an input a 
                % network adjacency matrix A, for a network that is undirected. If you
                % provide a network that is directed, this code is going to make it
                % undirected before continuing. Since link weights will not affect
                % component sizes, weighted and unweighted networks work equally well. You
                % may provide a "full" or a "sparse" matrix.
                %
                % OUTPUTS:
                % nComponents             INT - The number of components in the network.
                % sizes                 vector<INT> - a vector of component sizes, sorted, 
                %   descending.
                % members               cell<vector<INT>> a cell array of vectors, each
                %   entry of which is a membership list for that component, sorted, 
                %   descending by component size.
                
                % Number of nodes
                N = size(A,1);
                % Remove diagonals
                A(1:N+1:end) = 0;
                % make symmetric, just in case it isn't
                A=A+A';
                % Have we visited a particular node yet?
                isDiscovered = zeros(N,1);
                % Empty members cell
                members = {};
                % check every node
                for n=1:N
                    if ~isDiscovered(n)
                        % started a new group so add it to members
                        members{end+1} = n;
                        % account for discovering n
                        isDiscovered(n) = 1;
                        % set the ptr to 1
                        ptr = 1;
                        while (ptr <= length(members{end}))
                            % find neighbors
                            nbrs = find(A(:,members{end}(ptr)));
                            % here are the neighbors that are undiscovered
                            newNbrs = nbrs(isDiscovered(nbrs)==0);
                            % we can now mark them as discovered
                            isDiscovered(newNbrs) = 1;
                            % add them to member list
                            members{end}(end+1:end+length(newNbrs)) = newNbrs;
                            % increment ptr so we check the next member of this component
                            ptr = ptr+1;
                        end
                    end
                end
                % number of components
                nComponents = length(members);
                for n=1:nComponents
                    % compute sizes of components
                    sizes(n) = length(members{n});
                end

                [sizes,idx] = sort(sizes,'descend');
                members = members(idx);

            end
   
            
            function [d, vals] = infected_partner_ratio(state, adj)
           %% INFECTED_PARTNER_RATIO 
            % 
            %     state = ibm.indiv_state;
            %     adj = ibm.adj;

                s = any(state,2);

                d_seq = full(sum(adj,2));
                max_deg = max(d_seq);

                % only interested in row values for infected individuals
                adj(~s,:) = 0;

                % count only those with infected neighbours
                adj(:,~s) = 0;

                % col1: # infected partners,  col2: total # partners
                vals = [full(sum(adj(s,:),2)), d_seq(s)];

                d = zeros(1, max_deg);
                for i = 1:max_deg
                    d(i) = mean(vals(vals(:,2) == i,1));       
                end

                mean_ratio = mean(vals(:,1)./vals(:,2));


            end

            function [fig_h, plot_names] = plot_prev(data, day_range, DIR_NAME)
              %% Prevalence of each strain as percentage of population, over time
                % (static function can be called outside of any simulation)
                            
                N = data.N;
                days = day_range(1):day_range(2);
                plot_names = {};

                PREV_MAX = 20;
                
                FONT_SZ = 16;

                fig_h(1) = figure('name','Prevalence','color','w','position',[ 147 80 740 514]);
                subaxis(1,1,1, 'margin',0.15);
                set(gca,'fontsize',FONT_SZ,'linewidth',1)
                hold on;
                box on;
                
                if isfield(data,'prevalence_std')
                    
                    H(1) = shadedErrorBar(days,100*data.prev_either(day_range(1)+1:day_range(2)+1)./N,...
                        100*data.prev_either_std(day_range(1)+1:day_range(2)+1)./N,{'k-','linewidth',2},1);
                    
                    H(2) = shadedErrorBar(days,100*data.prevalence(day_range(1)+1:day_range(2)+1,1)./N,...
                        100*data.prevalence_std(day_range(1)+1:day_range(2)+1,1)./N,{'b-','linewidth',1},1);
                    
%                     H(2) = shadedErrorBar(days,100*data.prev_nonAMRonly(day_range(1)+1:day_range(2)+1,1)./N,...
%                         100*data.prev_nonAMRonly_std(day_range(1)+1:day_range(2)+1,1)./N,{'b-','linewidth',1},1);
                    
                    H(3) = shadedErrorBar(days,100*data.prevalence(day_range(1)+1:day_range(2)+1,2)./N,...
                        100*data.prevalence_std(day_range(1)+1:day_range(2)+1,2)./N,{'r-','linewidth',1},1);
                    
                    H(4) = shadedErrorBar(days,100*data.prev_both(day_range(1)+1:day_range(2)+1)./N,...
                        100*data.prev_both_std(day_range(1)+1:day_range(2)+1)./N,{'m--','linewidth',1},1);
                    
                    leg = legend([H(1).mainLine, H(2).mainLine, H(3).mainLine, H(4).mainLine],'either strain',...
                        'non-AMR strain','AMR strain','Coinfected');
                    set(leg,'location','northeast');
                    
                else
                
                    h_prev_either = plot(days, 100*data.prev_either(day_range(1)+1:day_range(2)+1)./N,'k-','linewidth',2);
                    %h_prev_nonAMR = plot(days, 100*data.prevalence(day_range(1)+1:day_range(2)+1,1)./N,'b-');
                    h_prev_nonAMR = plot(days, 100*data.prev_nonAMRonly(day_range(1)+1:day_range(2)+1,1)./N,'b-');
                    h_prev_AMR = plot(days, 100*data.prevalence(day_range(1)+1:day_range(2)+1,2)./N,'r-');
                    h_prev_both = plot(days, 100*data.prev_both(day_range(1)+1:day_range(2)+1)./N,'m--','linewidth',1);
                    

                    % legend
%                     leg = legend([h_prev_either,h_prev_nonAMR,h_prev_AMR,h_prev_both],'either strain',...
%                         'non-AMR strain','AMR strain','Coinfected','location','northwest');
                    
                    leg = legend([h_prev_either,h_prev_nonAMR,h_prev_AMR,h_prev_both],'either strain',...
                        'non-AMR strain (only)','AMR strain','Coinfected');
                    set(leg,'location','northeast');
                    
                end    
                axis([day_range(1) day_range(2) 0 PREV_MAX]);
                %grid on;
                box on;
                ylabel('Strain prevalence (%)','fontsize',FONT_SZ);
                xlabel('Time (days)','fontsize',FONT_SZ);
                set(leg,'fontsize',14,'linewidth',1);    
                title('Prevalence','fontsize',FONT_SZ);
                
                plot_names{1} = 'prevalence';

                % SAVE PLOT
                if ~isempty(DIR_NAME)
                    % check output folders exist, if not then create
                    if exist('./outdata','dir')~=7;mkdir('./outdata');end
                    if exist(['./outdata/',DIR_NAME],'dir')~=7;mkdir(['./outdata/',DIR_NAME]);end
                    export_fig(['./outdata/',DIR_NAME,'/prevalence.pdf'],fig_h(1));
                    %if isfield(data,'prevalence_std')
                    export_fig(['./outdata/',DIR_NAME,'/prevalence.png'],'-r150',fig_h(1));
                end
                %end
                
            end
            
            function [fig_h, plot_names] = plot_incidence(data, day_range, ave, DIR_NAME)
              %% Incidence of each strain - block averaged  (per week if ave=7) 
                % (static function can be called outside of any simulation)
                 
                N = data.N;
                
                % override = weekly block average
                ave = 7;
                day_range = [0 size(data.incidence,1)-1];
                
                days = day_range(1):day_range(2);
                plot_names = {};

                % average incidence over 'ave' blocks 
                % - express as per 1/denom
                denom = 10000;
                incidence(:,1) = arrayfun(@(i) mean(data.incidence_either(i:i+ave-1)),1:ave:length(data.incidence_either)-ave+1)'.*(denom/N);
                incidence(:,2) = arrayfun(@(i) mean(data.incidence(i:i+ave-1,1)),1:ave:length(data.incidence(:,1))-ave+1)'.*(denom/N);
                incidence(:,3) = arrayfun(@(i) mean(data.incidence(i:i+ave-1,2)),1:ave:length(data.incidence(:,2))-ave+1)'.*(denom/N);
    
                
                fig_h(1) = figure('name','Incidence','color','w','position',[200   270   655   503]);
                set(gca,'fontsize',12,'linewidth',1)
                hold on;
                box on;
                h_inc_either = plot(1:length(days)/ave,incidence(:,1),'k-','linewidth',2);
                h_inc_nonAMR = plot(1:length(days)/ave,incidence(:,2),'b-');
                h_inc_AMR = plot(1:length(days)/ave,incidence(:,3),'r-');
                xlim([day_range(1) day_range(2)]/ave);
                grid on;
                ylabel(['New infections: per day, per ' , num2str(denom),' people']);
                xlabel('Week');

                % legend
                leg = legend([h_inc_either,h_inc_nonAMR,h_inc_AMR],'either strain',...
                    'non-AMR strain','AMR strain','location','northwest');
                set(leg,'fontsize',11);
                title('Incidence');
                
                plot_names{1} = 'incidence';

                % SAVE PLOT
                % check output folders exist, if not then create
                if ~isempty(DIR_NAME)
                    if exist('./outdata','dir')~=7;mkdir('./outdata');end
                    if exist(['./outdata/',DIR_NAME],'dir')~=7;mkdir(['./outdata/',DIR_NAME]);end
                    export_fig(['./outdata/',DIR_NAME,'/incidence.pdf'],fig_h(1));
                    export_fig(['./outdata/',DIR_NAME,'/incidence.png'],'-r150',fig_h(1));
                end
                                
            end
            
            function [fig_h, plot_names] = plot_attendance(data, day_range, ave, DIR_NAME)
              %% Incidence of each strain - block averaged  (per week if ave=7) 
                % (static function can be called outside of any simulation)
                 
                N = data.N;
                
                % override = weekly block average
                ave = 7;
                day_range = [0 size(data.n_attended_seeked,1)-1];
                
                days = day_range(1):day_range(2);
                plot_names = {};

                
                fig_h(1) = figure('name','Attendance stats.','color','w','position',[110         426        1734         536]);
                
                % attended who were voluntary treatment seeking
                
                % average ivalues over 'ave' blocks
                seek(:,1) = arrayfun(@(i) mean(data.n_attended_seeked(i:i+ave-1,1)),1:ave:length(data.n_attended_seeked(:,1))-ave+1)';
                seek(:,2) = arrayfun(@(i) mean(data.n_attended_seeked(i:i+ave-1,2)),1:ave:length(data.n_attended_seeked(:,2))-ave+1)';
                seek(:,3) = arrayfun(@(i) mean(data.n_attended_seeked(i:i+ave-1,3)),1:ave:length(data.n_attended_seeked(:,3))-ave+1)';
    
                subaxis(1,3,1,'margin',0.1,'spacinghoriz',0.07)
                set(gca,'fontsize',12,'linewidth',1)
                hold on;
                box on;
%                 h_seeked = plot(days,data.n_attended_seeked(days+1,1),'k-','linewidth',2);
%                 h_seeked_nonAMR = plot(days,data.n_attended_seeked(days+1,2),'b-','linewidth',2);
%                 h_seeked_AMR = plot(days,data.n_attended_seeked(days+1,3),'r-','linewidth',2);
                
                h_seeked = plot(1:length(days)/ave,seek(:,1),'k-','linewidth',2);
                h_seeked_nonAMR = plot(1:length(days)/ave,seek(:,2),'b-','linewidth',1);
                h_seeked_AMR = plot(1:length(days)/ave,seek(:,3),'r-','linewidth',1);
                    
                %xlim([day_range(1) day_range(2)]);
                xlim([day_range(1) day_range(2)]/ave);
                grid on;
                ylabel(['# individuals attending per day']);
                %xlabel('Day');
                xlabel('Week');

                % legend
                leg = legend([h_seeked,h_seeked_nonAMR,h_seeked_AMR],'Any state',...
                    'non-AMR infected','AMR infected','location','northwest');
                set(leg,'fontsize',11);
                title('Attended voluntarily (SYMPTOMATIC)');
                
                
                % attended who were recalled
                
                % average values over 'ave' blocks
                recalled(:,1) = arrayfun(@(i) mean(data.n_attended_recalled(i:i+ave-1,1)),1:ave:length(data.n_attended_recalled(:,1))-ave+1)';
                recalled(:,2) = arrayfun(@(i) mean(data.n_attended_recalled(i:i+ave-1,2)),1:ave:length(data.n_attended_recalled(:,2))-ave+1)';
                recalled(:,3) = arrayfun(@(i) mean(data.n_attended_recalled(i:i+ave-1,3)),1:ave:length(data.n_attended_recalled(:,3))-ave+1)';
        
                
                subaxis(1,3,2)
                set(gca,'fontsize',12,'linewidth',1)
                hold on;
                box on;
%                 h_recalled = plot(days,data.n_attended_recalled(days+1,1),'k-','linewidth',2);
%                 h_recalled_nonAMR = plot(days,data.n_attended_recalled(days+1,2),'b-','linewidth',2);
%                 h_recalled_AMR = plot(days,data.n_attended_recalled(days+1,3),'r-','linewidth',2);

                h_recalled = plot(1:length(days)/ave,recalled(:,1),'k-','linewidth',2);
                h_recalled_nonAMR = plot(1:length(days)/ave,recalled(:,2),'b-','linewidth',1);
                h_recalled_AMR = plot(1:length(days)/ave,recalled(:,3),'r-','linewidth',1);

                %xlim([day_range(1) day_range(2)]);
                xlim([day_range(1) day_range(2)]/ave);
                
                grid on;
                %ylabel(['# individuals attending due to symptoms']);
                %xlabel('Day');
                xlabel('Week');
                title('Attended via RECALL');
                
                % attended who were traced
                
                % average values over 'ave' blocks
                traced(:,1) = arrayfun(@(i) mean(data.n_attended_traced(i:i+ave-1,1)),1:ave:length(data.n_attended_traced(:,1))-ave+1)';
                traced(:,2) = arrayfun(@(i) mean(data.n_attended_traced(i:i+ave-1,2)),1:ave:length(data.n_attended_traced(:,2))-ave+1)';
                traced(:,3) = arrayfun(@(i) mean(data.n_attended_traced(i:i+ave-1,3)),1:ave:length(data.n_attended_traced(:,3))-ave+1)';

                subaxis(1,3,3)
                set(gca,'fontsize',12,'linewidth',1)
                hold on;
                box on;
%                 h_traced = plot(days,data.n_attended_traced(days+1,1),'k-','linewidth',2);
%                 h_traced_nonAMR = plot(days,data.n_attended_traced(days+1,2),'b-','linewidth',2);
%                 h_traced_AMR = plot(days,data.n_attended_traced(days+1,3),'r-','linewidth',2);

                h_traced = plot(1:length(days)/ave,traced(:,1),'k-','linewidth',2);
                h_traced_nonAMR = plot(1:length(days)/ave,traced(:,2),'b-','linewidth',1);
                h_traced_AMR = plot(1:length(days)/ave,traced(:,3),'r-','linewidth',1);

                %xlim([day_range(1) day_range(2)]);
                xlim([day_range(1) day_range(2)]/ave);

                grid on;
                %ylabel(['# individuals attending due to symptoms']);
                %xlabel('Day');
                xlabel('Week');
                title('Attended via TRACE');
                
                
                plot_names{1} = 'attendance_stats';

                % SAVE PLOT
                if ~isempty(DIR_NAME)
                    % check output folders exist, if not then create
                    if exist('./outdata','dir')~=7;mkdir('./outdata');end
                    if exist(['./outdata/',DIR_NAME],'dir')~=7;mkdir(['./outdata/',DIR_NAME]);end
                    export_fig(['./outdata/',DIR_NAME,'/attendance_stats.pdf'],fig_h(1));
                    export_fig(['./outdata/',DIR_NAME,'/attendance_stats.png'],'-r150',fig_h(1));
                end
                
            end
            
            
            function [fig_h, plot_names] = plot_doses(data, ave, DIR_NAME)
              %% Plot number of doses of each drug and their effectiveness
                % (static function can be called outside of any simulation)
                   
                N = data.N;
                
                % override input
                ave = 30;
               
                 % y-axis maxima
                Y_MAX_CIPR = 300;  
                Y_MAX_CEFTA = 300; 
                
                plot_names = {};
                FONT_SZ = 20;
                
                % format data arrays appropriate for stacked bar chart
                cipr_data_full = [data.cipr_notinf, data.cipr_AMR, data.cipr_nonAMR];
                cefta_data_full = [data.cefta_notinf, data.cefta_nonAMR, data.cefta_AMR];
                rx_total_full = [data.cipr, data.cefta];
                cipr_data_ave = [];
                cefta_data_ave = [];
                rx_total_ave = [];
                
                % accumulate over number of days given by 'ave' (i.e.
                % weekly doses)
                n_vals = fix(size(cipr_data_full,1)/ave);
                for i = 1:n_vals
                    cipr_data_ave(i,:) = sum(cipr_data_full((i-1)*ave+1:i*ave,:),1);
                    cefta_data_ave(i,:) = sum(cefta_data_full((i-1)*ave+1:i*ave,:),1);             
                    rx_total_ave(i,:) = sum(rx_total_full((i-1)*ave+1:i*ave,:),1);             
                end
                
               
                
                fig_h(1) = figure('name','Doses','color','w','position',[ 331         209        1218         557]);
                % Cipr data
                subaxis(1,2,1,'margin',0.2,'spacinghoriz',0.1);
                set(gca,'fontsize',FONT_SZ,'linewidth',1)
                hold on;
                box on;
                b_cipr = bar(cipr_data_ave,'stacked');
                set(b_cipr,'edgecolor','none','linewidth',0.5,'barwidth',1);
                cols = {'r','y','g'};
                for n=1:3;set(b_cipr(n),'facecolor',cols{n});end
                P = findobj(gca,'type','patch');
                for n=1:length(P)
                    set(P(n),'facealpha',0.8)
                end
                axis([0 n_vals+1 0 Y_MAX_CIPR]);
                grid on;
                ylabel('# doses / month');
                if ave == 7
                    xlabel('Weeks','fontsize',FONT_SZ);
                elseif ave == 30
                    xlabel('Months','fontsize',FONT_SZ);
                else
                    xlabel(['Days / ',num2str(ave)],'fontsize',FONT_SZ);
                end
                title('Ciprofloxacin / Azithromycin','fontsize',20);
               
                % legend
                leg = legend('No infection (not required)','AMR infected (undertreated)',...
                    'Correct treatment (non-AMR)','location','northeast');
                set(leg,'fontsize',14,'linewidth',1);

                % CEFT/A data
                subaxis(1,2,2,'margin',0.2,'spacinghoriz',0.1);
                set(gca,'fontsize',FONT_SZ,'linewidth',1)
                hold on;
                box on;
                b_cefta = bar(cefta_data_ave,'stacked');
                set(b_cefta,'edgecolor','none','linewidth',0.5,'barwidth',1);
                cols = {'r','b','g'};
                for n=1:3;set(b_cefta(n),'facecolor',cols{n});end
                P = findobj(gca,'type','patch');
                for n=1:length(P)
                    set(P(n),'facealpha',0.8)
                end
                axis([0 n_vals+1 0 Y_MAX_CEFTA]);
                grid on;
                ylabel('# doses / month','fontsize',FONT_SZ);
                if ave == 7
                    xlabel('Weeks','fontsize',FONT_SZ);
                elseif ave == 30
                    xlabel('Months','fontsize',FONT_SZ);
                else
                    xlabel(['Days / ',num2str(ave)],'fontsize',FONT_SZ);
                end
                title('Ceftriaxone / Azithromycin','fontsize',20);
                
                % legend
                leg = legend('No infection (not required)','non-AMR infected (over treated)',...
                    'Correct treatment (AMR)','location','northeast');
                set(leg,'fontsize',14,'linewidth',1);
                        
                plot_names{1} = 'doses';
                
                % SAVE PLOT
                if ~isempty(DIR_NAME)
                    % check output folders exist, if not then create[h_fig_dose,plot_names] = AMR_IBM.plot_doses(data, N, ave, DIR_NAME);
                    if exist('./outdata','dir')~=7;mkdir('./outdata');end
                    if exist(['./outdata/',DIR_NAME],'dir')~=7;mkdir(['./outdata/',DIR_NAME]);end
                    export_fig(['./outdata/',DIR_NAME,'/dose_breakdown.pdf'],fig_h(1));
                    export_fig(['./outdata/',DIR_NAME,'/dose_breakdown.png'],'-r150',fig_h(1));
                end


                
              %% overall treatments (Rx)
                fig_h(2) = figure('name','#Rx','color','w','position',[564   199   713   586]);
                set(gca,'fontsize',16,'linewidth',1)
                hold on;
                box on;
                b_rx = bar(rx_total_ave,'stacked');
                set(b_rx,'edgecolor','none','linewidth',0.5,'barwidth',1);
                % legend
                leg = legend('Ciprofloxacin','Ceftriaxone','location','northeast');
                set(leg,'fontsize',11);
                cols = {'b','r'};
                for n=1:2;set(b_rx(n),'facecolor',cols{n});end                    
                P = findobj(gca,'type','patch');
                for n=1:length(P)
                    set(P(n),'facealpha',0.8)
                end
                axis([0 n_vals+1 0 max([Y_MAX_CIPR, Y_MAX_CEFTA])]);
                grid on;
                ylabel('#doses / month');
                if ave == 7
                    xlabel('Weeks');
                else
                    xlabel(['Days / ',num2str(ave)]);
                end
                title('# Rx','fontsize',16);
                
                % SAVE PLOT
                if ~isempty(DIR_NAME)
                    % check output folders exist, if not then create
                    if exist('./outdata','dir')~=7;mkdir('./outdata');end
                    if exist(['./outdata/',DIR_NAME],'dir')~=7;mkdir(['./outdata/',DIR_NAME]);end
                    export_fig(['./outdata/',DIR_NAME,'/Rx.pdf'],fig_h(2));
                    export_fig(['./outdata/',DIR_NAME,'/Rx.png'],'-r150',fig_h(2));
                end

                
                
                

                
%%               Treatment (doses), over time
%                 %self.fig_h(1) = figure('name','Treatment','color','w','position',[600 200 800 600]);
%                 subaxis(1,3,2);
%                 
%                 set(gca,'fontsize',16,'linewidth',1)
%                 hold on;
%                 box on;
%                 p_cipr = plot(days, smooth(data.cipr,7),'b-','linewidth',2);
%                 p_wcipr = plot(days, smooth(data.wasted_cipr,7),'b-');
%                 p_cefta = plot(days, smooth(data.cefta,7),'r-','linewidth',2);
%                 p_wcefta = plot(days, smooth(data.wasted_cefta,7),'r-');
%                 p_mistreated = plot(days, smooth(data.mistreated,7),'k-','linewidth',1);
%                 xlim([day_range(1) day_range(2)]);
%                 grid on;
%                 ylabel('Doses (per day)');
%                 xlabel('Time (days)');
%                 
%                 % legend
%                 leg = legend([p_cipr, p_cefta, p_mistreated],'Cipr (non-AMR) doses','Ceft/A (AMR) doses',...
%                      'AMR patient recalls','location','northwest');
%                 set(leg,'fontsize',11);
%                 
%                 plot_names{2} = 'treatments';
                
            end
            
            function [fig_h, plot_names] = plot_AMR_ratio(data, DIR_NAME)
              %% AMR strain ratio
              
                % compute the prevalence ratio of AMR to either (totalprev.)
                strain_ratio = data.prevalence(:,2)./data.prev_either;
                
                fig_h = figure('name','AMR strain ratio','color','w','position',[300          48         560         420]);
                set(gca,'fontsize',16)
                plot(strain_ratio*100,'k-');
                xlabel('Day');
                ylabel('% resistant isolates');
                
                plot_names = {'AMR_ratio'};
                
                % SAVE PLOT
                if ~isempty(DIR_NAME)
                    % check output folders exist, if not then create
                    if exist('./outdata','dir')~=7;mkdir('./outdata');end
                    if exist(['./outdata/',DIR_NAME],'dir')~=7;mkdir(['./outdata/',DIR_NAME]);end
                    export_fig(['./outdata/',DIR_NAME,'/AMR_ratio.pdf'],fig_h(1));
                    export_fig(['./outdata/',DIR_NAME,'/AMR_ratio.png'],'-r150',fig_h(1));
                end
    
            end
            
            function [fig_h, plot_names] = plot_diagnoses(data, DIR_NAME)
                
                % positive diagnoses (equal to number of doses given to
                % infecteds per day?)
                test_positive = (data.cefta - data.cefta_notinf) + ...
                                    (data.cipr - data.cipr_notinf);
                
                fig_h = figure('name','Diagnoses','color','w','position',[400          48         560         420]);
                set(gca,'fontsize',16);
                box on;
                grid on;
                hold on;
                h_incnew = plot(smooth(10000*data.incidence_new./data.N,7),'k-');
                h_tpos = plot(smooth(10000*test_positive./data.N,7),'r-');
                xlabel('Day');
                ylabel('# individuals / 10,000 / day');
                
                leg = legend([h_incnew,h_tpos],'Newly infected','Tested positive (either strain)');
                set(leg,'fontsize',12);
                
                plot_names = {'diagnoses'};
                
                
                % SAVE PLOT
                if ~isempty(DIR_NAME)
                    % check output folders exist, if not then create
                    if exist('./outdata','dir')~=7;mkdir('./outdata');end
                    if exist(['./outdata/',DIR_NAME],'dir')~=7;mkdir(['./outdata/',DIR_NAME]);end
                    export_fig(['./outdata/',DIR_NAME,'/diagnoses.pdf'],fig_h(1));
                    export_fig(['./outdata/',DIR_NAME,'/diagnoses.png'],'-r150',fig_h(1));
                end
            end
            
            function[fig_h, plot_name] = plot_dose_bar(data,window,DIR_NAME)
                %% 
                N = data.N;
               
                fig_h(1) = figure;
                shadedErrorBar(0:length(data.prev_either)-1,100*data.prev_either./N,100*data.prev_either_std./N,{'b'},1);
                
                
                plot_name = {'dose_bar'};
            end
            
                
            function [hfig,  h ] = plplot(x, xmin, alpha)
            % PLPLOT visualizes a power-law distributional model with empirical data.
            %    Source: http://www.santafe.edu/~aaronc/powerlaws/
            % 
            %    PLPLOT(x, xmin, alpha) plots (on log axes) the data contained in x 
            %    and a power-law distribution of the form p(x) ~ x^-alpha for 
            %    x >= xmin. For additional customization, PLPLOT returns a pair of 
            %    handles, one to the empirical and one to the fitted data series. By 
            %    default, the empirical data is plotted as 'bo' and the fitted form is
            %    plotted as 'k--'. PLPLOT automatically detects whether x is composed 
            %    of real or integer values, and applies the appropriate plotting 
            %    method. For discrete data, if min(x) > 50, PLFIT uses the continuous 
            %    approximation, which is a reliable in this regime.
            %
            %    Example:
            %       xmin  = 5;
            %       alpha = 2.5;
            %       x = xmin.*(1-rand(10000,1)).^(-1/(alpha-1));
            %       h = plplot(x,xmin,alpha);
            %
            %    For more information, try 'type plplot'
            %
            %    See also PLFIT, PLVAR, PLPVA

            % Version 1.0   (2008 February)
            % Copyright (C) 2008-2011 Aaron Clauset (Santa Fe Institute)
            % Distributed under GPL 2.0
            % http://www.gnu.org/copyleft/gpl.html
            % PLFIT comes with ABSOLUTELY NO WARRANTY
            % 
            % No notes
            % 
                % prepare figure window
                hfig = figure('color','w','position',[200 200 800 600]);
                set(gca,'fontsize',12);

                % reshape input vector
                x = reshape(x,numel(x),1);
                % initialize storage for output handles
                h = zeros(2,1);

                % select method (discrete or continuous) for plotting
                if     isempty(setdiff(x,floor(x))), f_dattype = 'INTS';
                elseif isreal(x),    f_dattype = 'REAL';
                else                 f_dattype = 'UNKN';
                end;
                if strcmp(f_dattype,'INTS') && min(x) > 50,
                    f_dattype = 'REAL';
                end;

                % estimate xmin and alpha, accordingly
                switch f_dattype,

                    case 'REAL',
                        n = length(x);
                        c = [sort(x) (n:-1:1)'./n];
                        q = sort(x(x>=xmin));
                        cf = [q (q./xmin).^(1-alpha)];
                        cf(:,2) = cf(:,2) .* c(find(c(:,1)>=xmin,1,'first'),2);

                        figure(hfig);
                        h(1) = loglog(c(:,1),c(:,2),'bo','MarkerSize',8,'MarkerFaceColor',[1 1 1]); hold on;
                        h(2) = loglog(cf(:,1),cf(:,2),'k--','LineWidth',2); hold off;
                        xr  = [10.^floor(log10(min(x))) 10.^ceil(log10(max(x)))];
                        xrt = (round(log10(xr(1))):2:round(log10(xr(2))));
                        if length(xrt)<4, xrt = (round(log10(xr(1))):1:round(log10(xr(2)))); end;
                        yr  = [10.^floor(log10(1/n)) 1];
                        yrt = (round(log10(yr(1))):2:round(log10(yr(2))));
                        if length(yrt)<4, yrt = (round(log10(yr(1))):1:round(log10(yr(2)))); end;
                        set(gca,'XLim',xr,'XTick',10.^xrt);
                        set(gca,'YLim',yr,'YTick',10.^yrt,'FontSize',16);
                        ylabel('Pr(X \geq x)','FontSize',16);
                        xlabel('x','FontSize',16)

                    case 'INTS',
                        n = length(x);        
                        q = unique(x);
                        c = hist(x,q)'./n;
                        c = [[q; q(end)+1] 1-[0; cumsum(c)]]; c(c(:,2)<10^-10,:) = [];
                        cf = ((xmin:q(end))'.^-alpha)./(zeta(alpha) - sum((1:xmin-1).^-alpha));
                        cf = [(xmin:q(end)+1)' 1-[0; cumsum(cf)]];
                        cf(:,2) = cf(:,2) .* c(c(:,1)==xmin,2);

                        figure(hfig);
                        h(1) = loglog(c(:,1),c(:,2),'bo','MarkerSize',8,'MarkerFaceColor',[1 1 1]); hold on;
                        h(2) = loglog(cf(:,1),cf(:,2),'k--','LineWidth',2); hold off;
                        xr  = [10.^floor(log10(min(x))) 10.^ceil(log10(max(x)))];
                        xrt = (round(log10(xr(1))):2:round(log10(xr(2))));
                        if length(xrt)<4, xrt = (round(log10(xr(1))):1:round(log10(xr(2)))); end;
                        yr  = [10.^floor(log10(1/n)) 1];
                        yrt = (round(log10(yr(1))):2:round(log10(yr(2))));
                        if length(yrt)<4, yrt = (round(log10(yr(1))):1:round(log10(yr(2)))); end;
                        set(gca,'XLim',xr,'XTick',10.^xrt);
                        set(gca,'YLim',yr,'YTick',10.^yrt,'FontSize',16);
                        ylabel('Pr(X \geq x)','FontSize',16);
                        xlabel('x','FontSize',16)

                    otherwise,
                        fprintf('(PLPLOT) Error: x must contain only reals or only integers.\n');
                        h = [];
                        return;
                end
            end
            
            function [alpha, xmin, L]=plfit(x, varargin)
            % PLFIT fits a power-law distributional model to data.
            %    Source: http://www.santafe.edu/~aaronc/powerlaws/
            % 
            %    PLFIT(x) estimates x_min and alpha according to the goodness-of-fit
            %    based method described in Clauset, Shalizi, Newman (2007). x is a 
            %    vector of observations of some quantity to which we wish to fit the 
            %    power-law distribution p(x) ~ x^-alpha for x >= xmin.
            %    PLFIT automatically detects whether x is composed of real or integer
            %    values, and applies the appropriate method. For discrete data, if
            %    min(x) > 1000, PLFIT uses the continuous approximation, which is 
            %    a reliable in this regime.
            %   
            %    The fitting procedure works as follows:
            %    1) For each possible choice of x_min, we estimate alpha via the 
            %       method of maximum likelihood, and calculate the Kolmogorov-Smirnov
            %       goodness-of-fit statistic D.
            %    2) We then select as our estimate of x_min, the value that gives the
            %       minimum value D over all values of x_min.
            %
            %    Note that this procedure gives no estimate of the uncertainty of the 
            %    fitted parameters, nor of the validity of the fit.
            %
            %    Example:
            %       x = (1-rand(10000,1)).^(-1/(2.5-1));
            %       [alpha, xmin, L] = plfit(x);
            %
            %    The output 'alpha' is the maximum likelihood estimate of the scaling
            %    exponent, 'xmin' is the estimate of the lower bound of the power-law
            %    behavior, and L is the log-likelihood of the data x>=xmin under the
            %    fitted power law.
            %    
            %    For more information, try 'type plfit'
            %
            %    See also PLVAR, PLPVA

            % Version 1.0    (2007 May)
            % Version 1.0.2  (2007 September)
            % Version 1.0.3  (2007 September)
            % Version 1.0.4  (2008 January)
            % Version 1.0.5  (2008 March)
            % Version 1.0.6  (2008 July)
            % Version 1.0.7  (2008 October)
            % Version 1.0.8  (2009 February)
            % Version 1.0.9  (2009 October)
            % Version 1.0.10 (2010 January)
            % Version 1.0.11 (2012 January)
            % Copyright (C) 2008-2012 Aaron Clauset (Santa Fe Institute)
            % Distributed under GPL 2.0
            % http://www.gnu.org/copyleft/gpl.html
            % PLFIT comes with ABSOLUTELY NO WARRANTY
            % 
            % Notes:
            % 
            % 1. In order to implement the integer-based methods in Matlab, the numeric
            %    maximization of the log-likelihood function was used. This requires
            %    that we specify the range of scaling parameters considered. We set
            %    this range to be [1.50 : 0.01 : 3.50] by default. This vector can be
            %    set by the user like so,
            %    
            %       a = plfit(x,'range',[1.001:0.001:5.001]);
            %    
            % 2. PLFIT can be told to limit the range of values considered as estimates
            %    for xmin in three ways. First, it can be instructed to sample these
            %    possible values like so,
            %    
            %       a = plfit(x,'sample',100);
            %    
            %    which uses 100 uniformly distributed values on the sorted list of
            %    unique values in the data set. Second, it can simply omit all
            %    candidates above a hard limit, like so
            %    
            %       a = plfit(x,'limit',3.4);
            %    
            %    Finally, it can be forced to use a fixed value, like so
            %    
            %       a = plfit(x,'xmin',3.4);
            %    
            %    In the case of discrete data, it rounds the limit to the nearest
            %    integer.
            % 
            % 3. When the input sample size is small (e.g., < 100), the continuous 
            %    estimator is slightly biased (toward larger values of alpha). To
            %    explicitly use an experimental finite-size correction, call PLFIT like
            %    so
            %    
            %       a = plfit(x,'finite');
            %    
            %    which does a small-size correction to alpha.
            %
            % 4. For continuous data, PLFIT can return erroneously large estimates of 
            %    alpha when xmin is so large that the number of obs x >= xmin is very 
            %    small. To prevent this, we can truncate the search over xmin values 
            %    before the finite-size bias becomes significant by calling PLFIT as
            %    
            %       a = plfit(x,'nosmall');
            %    
            %    which skips values xmin with finite size bias > 0.1.

                vec     = [];
                sample  = [];
                xminx   = [];
                limit   = [];
                finite  = false;
                nosmall = false;
                nowarn  = false;

                % parse command-line parameters; trap for bad input
                i=1; 
                while i<=length(varargin), 
                  argok = 1; 
                  if ischar(varargin{i}), 
                    switch varargin{i},
                        case 'range',        vec     = varargin{i+1}; i = i + 1;
                        case 'sample',       sample  = varargin{i+1}; i = i + 1;
                        case 'limit',        limit   = varargin{i+1}; i = i + 1;
                        case 'xmin',         xminx   = varargin{i+1}; i = i + 1;
                        case 'finite',       finite  = true;
                        case 'nowarn',       nowarn  = true;
                        case 'nosmall',      nosmall = true;
                        otherwise, argok=0; 
                    end
                  end
                  if ~argok, 
                    disp(['(PLFIT) Ignoring invalid argument #' num2str(i+1)]); 
                  end
                  i = i+1; 
                end
                if ~isempty(vec) && (~isvector(vec) || min(vec)<=1),
                    fprintf('(PLFIT) Error: ''range'' argument must contain a vector; using default.\n');
                    vec = [];
                end;
                if ~isempty(sample) && (~isscalar(sample) || sample<2),
                    fprintf('(PLFIT) Error: ''sample'' argument must be a positive integer > 1; using default.\n');
                    sample = [];
                end;
                if ~isempty(limit) && (~isscalar(limit) || limit<min(x)),
                    fprintf('(PLFIT) Error: ''limit'' argument must be a positive value >= 1; using default.\n');
                    limit = [];
                end;
                if ~isempty(xminx) && (~isscalar(xminx) || xminx>=max(x)),
                    fprintf('(PLFIT) Error: ''xmin'' argument must be a positive value < max(x); using default behavior.\n');
                    xminx = [];
                end;

                % reshape input vector
                x = reshape(x,numel(x),1);

                % select method (discrete or continuous) for fitting
                if     isempty(setdiff(x,floor(x))), f_dattype = 'INTS';
                elseif isreal(x),    f_dattype = 'REAL';
                else                 f_dattype = 'UNKN';
                end;
                if strcmp(f_dattype,'INTS') && min(x) > 1000 && length(x)>100,
                    f_dattype = 'REAL';
                end;

                % estimate xmin and alpha, accordingly
                switch f_dattype,

                    case 'REAL',
                        xmins = unique(x);
                        xmins = xmins(1:end-1);
                        if ~isempty(xminx),
                            xmins = xmins(find(xmins>=xminx,1,'first'));
                        end;
                        if ~isempty(limit),
                            xmins(xmins>limit) = [];
                        end;
                        if ~isempty(sample),
                            xmins = xmins(unique(round(linspace(1,length(xmins),sample))));
                        end;
                        dat   = zeros(size(xmins));
                        z     = sort(x);
                        for xm=1:length(xmins)
                            xmin = xmins(xm);
                            z    = z(z>=xmin); 
                            n    = length(z);
                            % estimate alpha using direct MLE
                            a    = n ./ sum( log(z./xmin) );
                            if nosmall,
                                if (a-1)/sqrt(n) > 0.1
                                    dat(xm:end) = [];
                                    xm = length(xmins)+1;
                                    break;
                                end;
                            end;
                            % compute KS statistic
                            cx   = (0:n-1)'./n;
                            cf   = 1-(xmin./z).^a;
                            dat(xm) = max( abs(cf-cx) );
                        end;
                        D     = min(dat);
                        xmin  = xmins(find(dat<=D,1,'first'));
                        z     = x(x>=xmin);
                        n     = length(z); 
                        alpha = 1 + n ./ sum( log(z./xmin) );
                        if finite, alpha = alpha*(n-1)/n+1/n; end; % finite-size correction
                        if n < 50 && ~finite && ~nowarn,
                            fprintf('(PLFIT) Warning: finite-size bias may be present.\n');
                        end;
                        L = n*log((alpha-1)/xmin) - alpha.*sum(log(z./xmin));

                    case 'INTS',

                        if isempty(vec),
                            vec  = (1.50:0.01:3.50);    % covers range of most practical 
                        end;                            % scaling parameters
                        zvec = zeta(vec);

                        xmins = unique(x);
                        xmins = xmins(1:end-1);
                        if ~isempty(xminx),
                            xmins = xmins(find(xmins>=xminx,1,'first'));
                        end;
                        if ~isempty(limit),
                            limit = round(limit);
                            xmins(xmins>limit) = [];
                        end;
                        if ~isempty(sample),
                            xmins = xmins(unique(round(linspace(1,length(xmins),sample))));
                        end;
                        if isempty(xmins)
                            fprintf('(PLFIT) Error: x must contain at least two unique values.\n');
                            alpha = NaN; xmin = x(1); D = NaN;
                            return;
                        end;
                        xmax   = max(x);
                        dat    = zeros(length(xmins),2);
                        z      = x;
                        fcatch = 0;

                        for xm=1:length(xmins)
                            xmin = xmins(xm);
                            z    = z(z>=xmin);
                            n    = length(z);
                            % estimate alpha via direct maximization of likelihood function
                            if fcatch==0
                                try
                                    % vectorized version of numerical calculation
                                    zdiff = sum( repmat((1:xmin-1)',1,length(vec)).^-repmat(vec,xmin-1,1) ,1);
                                    L = -vec.*sum(log(z)) - n.*log(zvec - zdiff);
                                catch
                                    % catch: force loop to default to iterative version for
                                    % remainder of the search
                                    fcatch = 1;
                                end;
                            end;
                            if fcatch==1
                                % force iterative calculation (more memory efficient, but 
                                % can be slower)
                                L       = -Inf*ones(size(vec));
                                slogz   = sum(log(z));
                                xminvec = (1:xmin-1);
                                for k=1:length(vec)
                                    L(k) = -vec(k)*slogz - n*log(zvec(k) - sum(xminvec.^-vec(k)));
                                end
                            end;
                            [Y,I] = max(L);
                            % compute KS statistic
                            fit = cumsum((((xmin:xmax).^-vec(I)))./ (zvec(I) - sum((1:xmin-1).^-vec(I))));
                            cdi = cumsum(hist(z,xmin:xmax)./n);
                            dat(xm,:) = [max(abs( fit - cdi )) vec(I)];
                        end
                        % select the index for the minimum value of D
                        [D,I] = min(dat(:,1));
                        xmin  = xmins(I);
                        z     = x(x>=xmin);
                        n     = length(z);
                        alpha = dat(I,2);
                        if finite, alpha = alpha*(n-1)/n+1/n; end; % finite-size correction
                        if n < 50 && ~finite && ~nowarn,
                            fprintf('(PLFIT) Warning: finite-size bias may be present.\n');
                        end;
                        L     = -alpha*sum(log(z)) - n*log(zvec(find(vec<=alpha,1,'last')) - sum((1:xmin-1).^-alpha));

                    otherwise,
                        fprintf('(PLFIT) Error: x must contain only reals or only integers.\n');
                        alpha = [];
                        xmin  = [];
                        L     = [];
                        return;
                end
            end
            
        end % static methods



end  % class