% script initialisation
clear;  clc;  close all;
figID = 0;

% Model and modality setting
model = 4;          % 1 = Hindmarsh-Rose model
                    % 2 = Epileptor model
                    % 3 = Jansen-Rit model
                    % 4 = simple autoactivating feedback loop
                    % 5 = Heat-shock model

modality = 4;       % 1 = gPC convergence plots
                    % 2 = comparison of different UQ methods
                    % 3 = nominal system phase portrait and deterministic parameter sweeping
                    % 4 = computation of the first two moments according to a specified surrogate model
                    % 5 = computation of execution times and RMSEs with the Galerkin approach
                    % 6 = computation of execution times and RMSEs with the collocation approach     

switch model 
    %% Hindmarsh-Rose model
    case 1
        % Initial conditions setting
        x_0 = [0,0,0];    % [0.3, 1, 0.2]
        
        % Model parameters setting
        a = 1.0;         % Parameter 1 (nominal 1.0, fixed)
        b = 2.525;       % Parameter 2 (nominal 3.0, plausible range [2, 5]) 
        c = 1.0;         % Parameter 3 (nominal 1.0, fixed)
        d = 5.0;         % Parameter 4 (nominal 5.0, fixed)
        I = 4.0;         % Parameter 5 (nominal 3.5, range [0.4, 5.6])
        r = 0.01;        % Parameter 6 (usual range [0.001, 0.035])
        s = 4.0;         % Parameter 7 (nominal 4.0, fixed)
        x_R = -8/5;      % Parameter 8 (nominal -8/5, fixed)
        
        parameters_vec = [a, b, c, d, I, r, s, x_R];
        
        % Setting of parameter that is uncertain around the nominal value
        parameter_of_interest = [2,5];    % parameter that will be uncertain, insert the parameter number
        uncertainty_delta = [0.05,0.4];   % nominal system value is at the center of this uncertainty interval range, order of uncertainty follows order of parameter_of_interest
        uncertainty = {parameter_of_interest,uncertainty_delta};
        
        % Parameter name assignment for plots
        parameter_name = parameter_name_assignment(parameter_of_interest,model);

        % Construction of structs for PoCET toolbox
        [states,parameters,inputs] = HR_struct_create(parameters_vec,x_0,uncertainty);


    %% Epileptor model
    case 2
        % Initial conditions setting
        x_0 = [0,-5,3,0,0,0];    % [0,-5,3,0,0,0]
        
        % Model parameters setting
        r_1 = -1.6;      % Parameter 1 (nominal -1.6, range [-3, 0])
        r_2 = 1;         % Parameter 2 (nominal 1, fixed) 
        I_1 = 3.1;       % Parameter 3 (nominal 3.1, range [1, 5])
        I_2 = 0.45;      % Parameter 4 (nominal 0.45, range [0, 0.5])
        tau_0 = 2857;    % Parameter 5 (nominal 2857, range [2000, 3000])
        m = 0;           % Parameter 6 (nominal 0, range [-16, 1.5])
        tau_2 = 10;      % Parameter 7 (nominal 10, fixed)
        gamma = 0.01;    % Parameter 8 (nominal 1/100, fixed)
        
        parameters_vec = [r_1, r_2, I_1, I_2, tau_0, m, tau_2, gamma];
        
        % Setting of parameter that is uncertain around the nominal value
        parameter_of_interest = [3];    % parameter that will be uncertain, insert the parameter number
        uncertainty_delta = [0.5];      % nominal system value is at the center of this uncertainty interval range, order of uncertainty follows order of parameter_of_interest
        uncertainty = {parameter_of_interest,uncertainty_delta};
        
        % Parameter name assignment for plots
        parameter_name = parameter_name_assignment(parameter_of_interest,model);
        
        % Construction of structs for PoCET toolbox
        [states,parameters,inputs] = Epileptor_struct_create(parameters_vec,x_0,uncertainty);


    %% Jansen-Rit model
    case 3
        % Initial conditions setting
        x_0 = [0,0,0,0,0,0];    
        
        % Model parameters setting
        C = 135;         % Parameter 1 (nominal 135, fixed)
        A = 3.25;        % Parameter 2 (nominal 3.25, fixed) 
        B = 22;          % Parameter 3 (nominal 22, fixed)
        a = 100;         % Parameter 4 (nominal 100, fixed)
        b = 50;          % Parameter 5 (nominal 50, fixed)
        V_0 = 6;         % Parameter 6 (nominal 6, fixed)
        nu_max = 5;      % Parameter 7 (nominal 5, fixed)
        r = 0.56;        % Parameter 8 (nominal 0.56, fixed)
        p = 150;         % Parameter 9 (range [0,400])

        parameters_vec = [C, A, B, a, b, V_0, nu_max, r, p];
        
        % Setting of parameter that is uncertain around the nominal value
        parameter_of_interest = [9];    % parameter that will be uncertain, insert the parameter number
        uncertainty_delta = [50];       % nominal system value is at the center of this uncertainty interval range, order of uncertainty follows order of parameter_of_interest
        uncertainty = {parameter_of_interest,uncertainty_delta};
        
        % Parameter name assignment for plots
        parameter_name = parameter_name_assignment(parameter_of_interest,model);
        
        % Construction of structs for PoCET toolbox
        [states,parameters,inputs] = JR_struct_create(parameters_vec,x_0,uncertainty);

    %% Simple autoactivating feedback loop
    case 4
        % Initial conditions setting
        x_0 = 0;
        
        % Model parameters setting
        n = 2.0;         % Parameter 1 (greater than 1)
        K = 0.15;        % Parameter 2 (range [0, 0.3]) 
        a = 2;         % Parameter 3 (range [1, 3])
        
        parameters_vec = [n, K, a];
        
        % Setting of parameter that is uncertain around the nominal value
        parameter_of_interest = [3];    % parameter that will be uncertain, insert the parameter number
        uncertainty_delta = [0.25];      % nominal system value is at the center of this uncertainty interval range, order of uncertainty follows order of parameter_of_interest
        uncertainty = {parameter_of_interest,uncertainty_delta};
        
        % Parameter name assignment for plots
        parameter_name = parameter_name_assignment(parameter_of_interest,model);

        % Construction of structs for PoCET toolbox
        [states,parameters,inputs] = autoactivating_feedback_struct_create(parameters_vec,x_0,uncertainty);


    %% Heat-shock model
    case 5
        
        % Initial conditions setting
        x_0 = [100000,1,0.1,0.05,10.5,1,0.0012,0.0002,0.0008,0.0036,0.0036,1];    
        
        % Model parameters setting
        k_p = 11.49;            % Parameter 1 (nominal 11.49)
        k_pp = 111.9;           % Parameter 2 (nominal 111.9) 
        k_s = 100.7;            % Parameter 3 (nominal 100.7)
        k_sp = 533.1;           % Parameter 4 (nominal 533.1)
        k_F = 0.9640;           % Parameter 5 (nominal 0.9640)
        k_Fp = 0.9557;          % Parameter 6 (nominal 0.9557)
        k_FG = 0.005574;        % Parameter 7 (nominal 0.005574)
        k_FGp = 0.08371;        % Parameter 8 (nominal 0.08371)
        k_FsG = 0.9131;         % Parameter 9 (nominal 0.9131)
        k_FsGp = 0.4136;        % Parameter 10 (nominal 0.4136)
        k_Fs = 0.01192;         % Parameter 11 (nominal 0.01192)
        k_Fsp = 0.01064;        % Parameter 12 (nominal 0.01064)
        k_piRF = 18.16;         % Parameter 13 (nominal 18.16)
        k_piRH = 4.193;         % Parameter 14 (nominal 4.193)
        k_piHP = 0.5142;        % Parameter 15 (nominal 0.5142)
        k_piF = 0.02112;        % Parameter 16 (nominal 0.02112)
        d_F = 0.0008728;        % Parameter 17 (nominal 0.0008728)
        d_HP = 0.0009384;       % Parameter 18 (nominal 0.0009384)
        d_RF = 0.001719;        % Parameter 19 (nominal 0.001719)
        d_RP = 0.001017;        % Parameter 20 (nominal 0.001017)
        Ea = 174440;            % Parameter 21 (nominal 174440, range [160000, 830000]
        T_C = 36;               % Parameter 22 (nominal 36, range [20, 40]
        
        parameters_vec = [k_p,k_pp,k_s,k_sp,k_F,k_Fp,k_FG,k_FGp,k_FsG,k_FsGp,k_Fs,k_Fsp,k_piRF,k_piRH,k_piHP,k_piF,d_F,d_HP,d_RF,d_RP,Ea,T_C];
        
        % Setting of parameter that is uncertain around the nominal value
        parameter_of_interest = [22];    % parameter that will be uncertain, insert the parameter number
        uncertainty_delta = [20];        % nominal system value is at the center of this uncertainty interval range, order of uncertainty follows order of parameter_of_interest
        uncertainty = {parameter_of_interest,uncertainty_delta};
        
        % Parameter name assignment for plots
        parameter_name = parameter_name_assignment(parameter_of_interest,model);
        
        % Construction of structs for PoCET toolbox
        [states,parameters,inputs] = HS_struct_create(parameters_vec,x_0,uncertainty);

end

%% Simulation settings
t_initial = 0;
t_final = 500;
delta_t = 0.01;
integrator = 'ode45';
n_sweeping = 5;           % number of plots for deterministic_parameter_sweeping()
threshold_abs = 0.45;     % absolute threshold to flag disagreement between the ground truth and the UQ methods
threshold_rel = 0.02;     % relative threshold in percentage to flag disagreement between the ground truth and the UQ methods

% Statistics setting
pce_order = 7;
mc_samples = 500;         % for ground truth simulations select 1e5 samples
colloc_samples = 400;

% Options for the functions to select the desired surrogate model
galerkin = 1;
collocation = 2;
monte_carlo = 3;
galerkin_vs_mc = 1;
collocation_vs_mc = 2;
galerkin_vs_collocation = 3;

% Overall simulation options and uncertainty settings
simulation_opts = {t_initial,t_final,delta_t,integrator,mc_samples,colloc_samples, ...
                    pce_order,n_sweeping,threshold_abs,threshold_rel};

%% Computation of the desired performance

switch modality

    case 1
        % Loading of the ground truth (delta_t = 0.01, 1e5 Monte Carlo samples, t_final = 600)
        if (parameter_of_interest(1) == 2 && b == 2.25 && (b+uncertainty_delta(1)/2) == 2.5 && (b-uncertainty_delta(1)/2) == 2)
            ground_truth = load('ground truth b in 2_0 2_5');
        elseif (parameter_of_interest(1) == 2 && b == 2.75 && (b+uncertainty_delta(1)/2) == 3 && (b-uncertainty_delta(1)/2) == 2.5)
            ground_truth = load('ground truth b in 2_5 3_0');
        end

        max_order = 22;     % maximum expansion order to reach

        figID = gPC_convergence(states,parameters,inputs,ground_truth, max_order,collocation,simulation_opts,figID); 

    case 2
        % Loading of the ground truth (delta_t = 0.01, 1e5 Monte Carlo samples, t_final = 600)
        if (parameter_of_interest(1) == 2 && b == 2.25 && (b+uncertainty_delta(1)/2) == 2.5 && (b-uncertainty_delta(1)/2) == 2)
            ground_truth = load('ground truth b in 2_0 2_5');
        elseif (parameter_of_interest(1) == 2 && b == 2.75 && (b+uncertainty_delta(1)/2) == 3 && (b-uncertainty_delta(1)/2) == 2.5)
            ground_truth = load('ground truth b in 2_5 3_0');
        end

        [figID,times,~,~] = UQ_comparison(states,parameters,inputs,ground_truth,galerkin_vs_collocation,simulation_opts,figID);




    case 3
        % Nominal system simulation
        simoptions.tspan = [t_initial, t_final];
        simoptions.dt = delta_t;
        simoptions.setup = odeset;
        simoptions.solver = integrator; 
        t_span = [t_initial, t_final];
        
        if model == 1
            [t,Y] = ode45(@(t,Y)Hindmarsh_Rose(t,Y,parameters_vec),t_span,x_0,simoptions);
    
            figID = figID+1;
            figure(figID);
            plot(t,Y(:,1),'LineWidth',1.5);
            hold on;
            plot(t,x_R*ones(1,length(Y(:,1))),'LineWidth',1.5);  % resting potential line
            % hold on;
            % stem(t,2.5*peaks);
            % hold on;
            % plot(t,direct_current,'r','LineWidth',1.5);   % do not consider for the moment
            title('Membrane potential');
            xlabel('time','FontName','Arial','FontSize',14);  
            ylabel('x','FontName','Arial','FontSize',14);
            legend({'x','x_R','local max'});
            
            figID = figID+1;
            figure(figID);
            plot3(Y(:,1),Y(:,2),Y(:,3),'LineWidth',1.5);
            title('Phase portrait');
            xlabel('x','FontName','Arial','FontSize',14);  
            ylabel('y','FontName','Arial','FontSize',14);  
            zlabel('z','FontName','Arial','FontSize',14);
            grid on;
    
            % Deterministic parameter sweeping
            figID = deterministic_parameter_sweeping(parameters_vec,x_0,uncertainty,simulation_opts,figID,model);
            
            cut_off_indx = find(t > 600,1);
            [activity1,figID] = HR_neuronal_activity(Y(cut_off_indx:end,1),t(cut_off_indx:end),x_R,figID);

        elseif model == 2
            [t,Y] = ode45(@(t,Y)Epileptor(t,Y,parameters_vec),t_span,x_0,simoptions);

            figID = figID+1;
            figure(figID);
            plot(t,Y(:,1),'LineWidth',1.5);
            xlabel('time','FontName','Arial','FontSize',14);  
            ylabel('x_1','FontName','Arial','FontSize',14);

            figID = figID+1;
            figure(figID);
            plot(t,Y(:,4),'LineWidth',1.5);
            xlabel('time','FontName','Arial','FontSize',14);  
            ylabel('x_4','FontName','Arial','FontSize',14);

            figID = figID+1;
            figure(figID);
            plot(t,-Y(:,1)+Y(:,4),'LineWidth',1.5);
            xlabel('time','FontName','Arial','FontSize',14);  
            ylabel('-x_1+x_4','FontName','Arial','FontSize',14);
            
            figID = figID+1;
            figure(figID);
            plot3(-Y(:,1),Y(:,2),Y(:,3),'LineWidth',1.5);
            xlabel('-x_1','FontName','Arial','FontSize',14);  
            ylabel('x_2','FontName','Arial','FontSize',14);  
            zlabel('x_3','FontName','Arial','FontSize',14);
            grid on;

            % Deterministic parameter sweeping
            figID = deterministic_parameter_sweeping(parameters_vec,x_0,uncertainty,simulation_opts,figID,model);

       elseif model == 3
            [t,Y] = ode45(@(t,Y)Jansen_Rit(t,Y,parameters_vec),t_span,x_0,simoptions);

            figID = figID+1;
            figure(figID);
            plot(t,Y(:,2)-Y(:,3),'LineWidth',1.5);
            xlabel('time','FontName','Arial','FontSize',14);  
            ylabel('y_2-y_3','FontName','Arial','FontSize',14);

            % Deterministic parameter sweeping
            figID = deterministic_parameter_sweeping(parameters_vec,x_0,uncertainty,simulation_opts,figID,model);
            
        elseif model == 4
            [t,Y] = ode45(@(t,Y)autoactivating_feedback_loop(t,Y,parameters_vec),t_span,x_0,simoptions);

            figID = figID+1;
            figure(figID);
            plot(t,Y,'LineWidth',1.5);
            xlabel('time','FontName','Arial','FontSize',14);  
            ylabel('x','FontName','Arial','FontSize',14);

            % Deterministic parameter sweeping
            figID = deterministic_parameter_sweeping(parameters_vec,x_0,uncertainty,simulation_opts,figID,model);

        end








    case 4
        % Plot membrane potential mean of surrogate model
        surr_model = collocation;
        ground_truth = load('C) mean potential field');     % to be changed according to the regime analysed
        ground_truth = ground_truth.mean_potential;
        cut_off_indx = find(ground_truth.time > 600,1);
        
        if (model == 2 && surr_model == galerkin) ||  (model == 3 && surr_model == galerkin) ||  (model == 4 && surr_model == galerkin)
            disp('Galerkin approach not available for the selected model');
            return;
        end
        [results,simulation_time] = surrogate_model(states,parameters,inputs,surr_model,simulation_opts,model);
        
        if model == 1
            peaks1 = islocalmax(results.x.moments(1,:),'MinProminence',0.07);     % unfiltered data
            mean_potential_filt = sgolayfilt(results.x.moments(1,:),3,599);
            peaks2 = islocalmax(mean_potential_filt,'MinProminence',0.07);               % Savitzky-Golay filtering
    
            % Plot of unfiletered mean membrane potential
            figID = figID+1;
            figure(figID);
            plot(results.time,results.x.moments(1,:),'LineWidth',1.5);
            hold on;
            plot(results.time,parameters_vec(end)*ones(1,length(results.x.moments(1,:))),'LineWidth',1.5);  % resting potential line
            hold on;
            stem(results.time,2.5*peaks1);
    
            title_string = 'Mean membrane potential with ';
            for i=1:length(parameter_name)
                title_string = strcat([title_string,parameter_name(i),' \in ', ...
                    '[',num2str(parameters_vec(parameter_of_interest(i))-uncertainty_delta(i)/2),',',num2str(parameters_vec(parameter_of_interest(i))+uncertainty_delta(i)/2),']',', ']);
            end
            title_string = string(title_string(1:end-1));
            title(join(title_string));
    
            xlabel('time','FontName','Arial','FontSize',14);  
            ylabel('E[x]','FontName','Arial','FontSize',14);
    
            % Plot of filtered mean membrane potential
            figID = figID+1;
            figure(figID);
            plot(results.time,mean_potential_filt,'LineWidth',1.5);
            hold on;
            plot(results.time,parameters_vec(end)*ones(1,length(results.x.moments(1,:))),'LineWidth',1.5);  % resting potential line
            hold on;    
            stem(results.time,2.5*peaks2);
    
            title_string = 'Filtered mean membrane potential with ';
            for i=1:length(parameter_name)
                title_string = strcat([title_string,parameter_name(i),' \in ', ...
                    '[',num2str(parameters_vec(parameter_of_interest(i))-uncertainty_delta(i)/2),',',num2str(parameters_vec(parameter_of_interest(i))+uncertainty_delta(i)/2),']',', ']);
            end
            title_string = string(title_string(1:end-1));
            title(join(title_string));
            xlabel('time','FontName','Arial','FontSize',14);  
            ylabel('E[x]','FontName','Arial','FontSize',14);
    
            % Compute RMSE
            rmse_mean = rms(results.x.moments(1,cut_off_indx:end)-ground_truth.x.moments(1,cut_off_indx:end));
            rmse_variance = rms(results.x.moments(2,cut_off_indx:end)-ground_truth.x.moments(2,cut_off_indx:end));
    
            % Plot of histograms
            [activity2,figID] = HR_neuronal_activity(results.x.moments(1,cut_off_indx:end),results.time(cut_off_indx:end),x_R,figID);

        elseif model == 2
            
            % Plot x_1
            figID = figID+1;
            figure(figID);
            plot(results.time,results.x_1.moments(1,:),'LineWidth',1.5);
            title_string = 'Mean x_1 with ';
            for i=1:length(parameter_name)
                title_string = strcat([title_string,parameter_name(i),' \in ', ...
                    '[',num2str(parameters_vec(parameter_of_interest(i))-uncertainty_delta(i)/2),',',num2str(parameters_vec(parameter_of_interest(i))+uncertainty_delta(i)/2),']',', ']);
            end
            title_string = string(title_string(1:end-1));
            title(join(title_string));
    
            xlabel('time','FontName','Arial','FontSize',14);  
            ylabel('E[x_1]','FontName','Arial','FontSize',14);
            
            % Plot x_4
            figID = figID+1;
            figure(figID);
            plot(results.time,results.x_4.moments(1,:),'LineWidth',1.5);
            title_string = 'Mean x_4 with ';
            for i=1:length(parameter_name)
                title_string = strcat([title_string,parameter_name(i),' \in ', ...
                    '[',num2str(parameters_vec(parameter_of_interest(i))-uncertainty_delta(i)/2),',',num2str(parameters_vec(parameter_of_interest(i))+uncertainty_delta(i)/2),']',', ']);
            end
            title_string = string(title_string(1:end-1));
            title(join(title_string));
    
            xlabel('time','FontName','Arial','FontSize',14);  
            ylabel('E[x_4]','FontName','Arial','FontSize',14);
            
            % Plot -x_1+x_4
            figID = figID+1;
            figure(figID);
            plot(results.time,-results.x_1.moments(1,:)+results.x_4.moments(1,:),'LineWidth',1.5);
            title_string = 'Mean -x_1+x_4 with ';
            for i=1:length(parameter_name)
                title_string = strcat([title_string,parameter_name(i),' \in ', ...
                    '[',num2str(parameters_vec(parameter_of_interest(i))-uncertainty_delta(i)/2),',',num2str(parameters_vec(parameter_of_interest(i))+uncertainty_delta(i)/2),']',', ']);
            end
            title_string = string(title_string(1:end-1));
            title(join(title_string));
    
            xlabel('time','FontName','Arial','FontSize',14);  
            ylabel('E[-x_1+x_4]','FontName','Arial','FontSize',14);
    
        elseif model == 3
            % Plot y_2
            figID = figID+1;
            figure(figID);
            plot(results.time,results.y_2.moments(1,:),'LineWidth',1.5);
            title_string = 'Mean y_2 with ';
            for i=1:length(parameter_name)
                title_string = strcat([title_string,parameter_name(i),' \in ', ...
                    '[',num2str(parameters_vec(parameter_of_interest(i))-uncertainty_delta(i)/2),',',num2str(parameters_vec(parameter_of_interest(i))+uncertainty_delta(i)/2),']',', ']);
            end
            title_string = string(title_string(1:end-1));
            title(join(title_string));
    
            xlabel('time','FontName','Arial','FontSize',14);  
            ylabel('E[y_2]','FontName','Arial','FontSize',14);

            % Plot y_3
            figID = figID+1;
            figure(figID);
            plot(results.time,results.y_3.moments(1,:),'LineWidth',1.5);
            title_string = 'Mean y_3 with ';
            for i=1:length(parameter_name)
                title_string = strcat([title_string,parameter_name(i),' \in ', ...
                    '[',num2str(parameters_vec(parameter_of_interest(i))-uncertainty_delta(i)/2),',',num2str(parameters_vec(parameter_of_interest(i))+uncertainty_delta(i)/2),']',', ']);
            end
            title_string = string(title_string(1:end-1));
            title(join(title_string));
    
            xlabel('time','FontName','Arial','FontSize',14);  
            ylabel('E[y_3]','FontName','Arial','FontSize',14);

            % Plot y_2 - y_3
            figID = figID+1;
            figure(figID);
            plot(results.time,results.y_2.moments(1,:) - results.y_3.moments(1,:),'LineWidth',1.5);
            title_string = 'Mean y_2 - y_3 with ';
            for i=1:length(parameter_name)
                title_string = strcat([title_string,parameter_name(i),' \in ', ...
                    '[',num2str(parameters_vec(parameter_of_interest(i))-uncertainty_delta(i)/2),',',num2str(parameters_vec(parameter_of_interest(i))+uncertainty_delta(i)/2),']',', ']);
            end
            title_string = string(title_string(1:end-1));
            title(join(title_string));
    
            xlabel('time','FontName','Arial','FontSize',14);  
            ylabel('E[y_2 - y_3]','FontName','Arial','FontSize',14);

        elseif model == 4
            % Plot mean of x
            figID = figID+1;
            figure(figID);
            plot(results.time,results.x.moments(1,:),'LineWidth',1.5);
            title_string = 'Mean x with ';
            for i=1:length(parameter_name)
                title_string = strcat([title_string,parameter_name(i),' \in ', ...
                    '[',num2str(parameters_vec(parameter_of_interest(i))-uncertainty_delta(i)/2),',',num2str(parameters_vec(parameter_of_interest(i))+uncertainty_delta(i)/2),']',', ']);
            end
            title_string = string(title_string(1:end-1));
            title(join(title_string));
    
            xlabel('time','FontName','Arial','FontSize',14);  
            ylabel('E[x]','FontName','Arial','FontSize',14);

            % Plot variance of x
            figID = figID+1;
            figure(figID);
            plot(results.time,results.x.moments(2,:),'LineWidth',1.5);
            title_string = 'Variance x with ';
            for i=1:length(parameter_name)
                title_string = strcat([title_string,parameter_name(i),' \in ', ...
                    '[',num2str(parameters_vec(parameter_of_interest(i))-uncertainty_delta(i)/2),',',num2str(parameters_vec(parameter_of_interest(i))+uncertainty_delta(i)/2),']',', ']);
            end
            title_string = string(title_string(1:end-1));
            title(join(title_string));
    
            xlabel('time','FontName','Arial','FontSize',14);  
            ylabel('Var[x]','FontName','Arial','FontSize',14);
            
        end














    case 5
        max_expansion_order = 35;
        
        % Variables to be saved after completing simulation
        execution_time_galerkin = zeros(1,max_expansion_order);
        rmse_mean = zeros(1,max_expansion_order);
        rmse_variance = zeros(1,max_expansion_order);

        ground_truth = load('D) mean potential field');     % to be changed according to the regime analysed
        ground_truth = ground_truth.mean_potential;

        cut_off_indx = find(ground_truth.time > 600,1);

        for i=1:max_expansion_order
            simulation_opts{1,7} = i;
            [mean_potential,simulation_time_G] = HR_surrogate_model(states,parameters,galerkin,simulation_opts);

            execution_time_galerkin(i) = simulation_time_G;
            rmse_mean(i) = rms(mean_potential.x.moments(1,cut_off_indx:end)-ground_truth.x.moments(1,cut_off_indx:end));
            rmse_variance(i) = rms(mean_potential.x.moments(2,cut_off_indx:end)-ground_truth.x.moments(2,cut_off_indx:end));
        end
        
        % Plot of simulation time for Galerkin with increasing order of expansion
        X = 1:max_expansion_order;
        fit1 = fit(X',execution_time_galerkin','exp1');
        fit2 = fit(X',execution_time_galerkin','exp2');

        figID = figID+1;
        figure(figID);
        plot(execution_time_galerkin,'LineWidth',1.5);
        hold on;
        plot(fit2,X',execution_time_galerkin');
        set(gca, 'YScale','log');
        title('Execution time Galerkin approach');
        xlabel('Expansion order');
        ylabel('Logarithm of execution time [s]');
        
         % Plot of RMSEs for Galerkin with increasing order of expansion
        figID = figID+1;
        figure(figID);
        plot(rmse_mean,'LineWidth',1.5);
        title('RMSE mean');
        xlabel('Expansion order');
        ylabel('Mean RMSE');

        figID = figID+1;
        figure(figID);
        plot(rmse_variance,'LineWidth',1.5);
        title('RMSE variance');
        xlabel('Expansion order');
        ylabel('Variance RMSE');

    case 6

        max_expansion_order = 12;
        max_colloc_samples = 50;
        
        % variables to be saved
        execution_time_collocation = zeros(max_expansion_order,int32(max_colloc_samples/50));
        rmse_mean = zeros(max_expansion_order,int32(max_colloc_samples/50));
        rmse_variance = zeros(max_expansion_order,int32(max_colloc_samples/50));

        ground_truth = load('C) mean potential field');     % to be changed according to the regime analysed
        ground_truth = ground_truth.mean_potential;

        cut_off_indx = find(ground_truth.time > 600,1);

        for i=1:max_expansion_order
            for j=1:int32(max_colloc_samples/50)
                simulation_opts{1,7} = i;
                simulation_opts{1,6} = j*50;
                [mean_potential,simulation_time_C] = HR_surrogate_model(states,parameters,collocation,simulation_opts);
                
                % expansion order increases from the last row to the first
                % row, while number of samples increases along the column
                % number
                execution_time_collocation(end-i+1,j) = simulation_time_C; 
                rmse_mean(end-i+1,j) = rms(mean_potential.x.moments(1,cut_off_indx:end)-ground_truth.x.moments(1,cut_off_indx:end));
                rmse_variance(end-i+1,j) = rms(mean_potential.x.moments(2,cut_off_indx:end)-ground_truth.x.moments(2,cut_off_indx:end));

                clear simulation_time_C;
            end
        end
        
        % Plot of the heatmaps
        figID = figID+1;
        figure(figID);
        xvalues = 50:50:max_colloc_samples;
        yvalues = 1:max_expansion_order;
        yvalues = yvalues(end:-1:1);    % for sake of plot visualisation
        h1 = heatmap(xvalues,yvalues,execution_time_collocation);
        h1.Title = 'Execution time collocation approach';
        h1.XLabel = 'Collocation samples';
        h1.YLabel = 'Expansion order';

        figID = figID+1;
        figure(figID);
        h2 = heatmap(xvalues,yvalues,rmse_mean);
        h2.Title = 'RMSE mean collocation approach';
        h2.XLabel = 'Collocation samples';
        h2.YLabel = 'Expansion order';

        figID = figID+1;
        figure(figID);
        h3 = heatmap(xvalues,yvalues,rmse_variance);
        h3.Title = 'RMSE variance collocation approach';
        h3.XLabel = 'Collocation samples';
        h3.YLabel = 'Expansion order';


        % ------------------------------------------------------------------------------------------------------------------------------------
        % % Plot of neuronal activity comparisons
        % filt = 1;
        % unfilt = 0;
        % figID = HR_activity_on_interval(parameters_vec,x_0,uncertainty,mean_potential,surr_model,unfilt,simulation_opts,figID);
        % figID = HR_activity_on_interval(parameters_vec,x_0,uncertainty,mean_potential,surr_model,filt,simulation_opts,figID);

        % % Computation of performance indices
        % [var1] = HR_neuronal_activity(Y(:,1),t,parameters_vec(end));

        % % % Plot the inter-spike interval
        % % figID = figID+1;
        % % figure(figID);
        % % stairs(var1.ISI,'LineWidth',1.5);
        % % title(strcat(['Membrane potential with b = ',num2str(parameters_vec(2)),', I = ',num2str(parameters_vec(5))])); 
        % % ylabel('ISI','FontName','Arial','FontSize',14);
        % 
        % % Plot the spike-count
        % figID = figID+1;
        % figure(figID);
        % stem(var1.spike_count,'LineWidth',1.5);
        % title(strcat(['Membrane potential with b = ',num2str(parameters_vec(2)),', I = ',num2str(parameters_vec(5))])); 
        % xlabel('Neuronal activity cycle','FontName','Arial','FontSize',14)
        % ylabel('Spike Count','FontName','Arial','FontSize',14);

end


