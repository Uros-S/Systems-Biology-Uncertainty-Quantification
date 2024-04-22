function figID = HR_activity_on_interval(nom_parameters,x_0,uncertainty,mean_potential,surrogate,filt,simulation_opts,figID)
% This function plots the spike count, duty cycle and ISI of some
% deterministic membrane potentials in a certain parameter range. On top of
% those plots, the same indices are plot considering the mean membrane
% potential obtained with a specified surrogate model.

    % Setting of simulation parameters
    parameter_number = uncertainty{1,1};
    interval = uncertainty{1,2};
    resolution = simulation_opts{1,8};
    t_span = [simulation_opts{1,1}, simulation_opts{1,2}];
    simoptions.dt = simulation_opts{1,3};
    simoptions.setup = odeset;
    simoptions.solver = {1,4};
    x_R = nom_parameters(end);

    % Filtering
    if filt
        mean_potential.x.moments(1,:) = sgolayfilt(mean_potential.x.moments(1,:),3,599);
    end

    % Name assignments for x-axis of the plot
    if parameter_number == 1
        parameter_name = 'a';
    elseif parameter_number == 2
        parameter_name = 'b';
    elseif parameter_number == 3
        parameter_name = 'c';
    elseif parameter_number == 4
        parameter_name = 'd';
    elseif parameter_number == 5
        parameter_name = 'I';
    elseif parameter_number == 6
        parameter_name = 'r';
    elseif parameter_number == 7
        parameter_name = 's';
    else
        parameter_name = 'x_R';
    end

    % Name assignment for the surrogate model
    if surrogate == 1
        pce_order = simulation_opts{1,7};
        model = strcat(['Galerkin method with order ',num2str(pce_order),' expansion']);
    elseif surrogate == 2
        pce_order = simulation_opts{1,7};
        colloc_samples = simulation_opts{1,6};
        model = strcat(['Collocation method with ',num2str(colloc_samples),' samples and order ',num2str(pce_order),' expansion']);
    elseif surrogate == 3
        mc_samples = simulation_opts{1,5};
        model = strcat(['Monte Carlo method with ',num2str(mc_samples),' samples']);
    end

    % Creation of a vector of different values of a parameter of interest
    parameter = nom_parameters(parameter_number); 
    parameter_lower = parameter-interval/2; 
    parameter_upper = parameter+interval/2;
    parameter_vec = parameter_lower:resolution:parameter_upper;
    
    % Creation of a matrix of parameters in which only the parameter of interest changes
    parameter_mat = kron(nom_parameters',ones(1,length(parameter_vec)));
    parameter_mat(parameter_number,:) = parameter_vec;

    % Setting lambda function for plot of membrane potential evolution w.r.t. a changing parameter
    lightBLUE = [0.356862745098039,0.811764705882353,0.956862745098039];
    darkBLUE = [0.0196078431372549,0.0745098039215686,0.670588235294118];
    blueGRADIENTflexible = @(i,N) lightBLUE + (darkBLUE-lightBLUE)*((i-1)/(N-1));

    % Compute the neuronal activity of the mean potential
    neur_act_mean = HR_neuronal_activity(mean_potential.x.moments(1,:),mean_potential.time,x_R);
    
    %% Plot of the inter-spike interval
    flag = 1;
    
    % Sequence of simulations in the parameter range of interest
    for i=1:length(parameter_vec)

        if flag                 % Initialisation of figure for plot of membrane potential
            figID = figID+1;
            figure(figID);
            p3 = stairs(neur_act_mean.ISI,'r');
            p3.LineWidth = 2;
            hold on;
            flag = 0;
        end

        [t,Y] = ode45(@(t,Y)Hindmarsh_Rose(t,Y,parameter_mat(:,i)),t_span,x_0,simoptions);
        [neur_act] = HR_neuronal_activity(Y(:,1),t,x_R);

        if i == 1                                % for inserting in the plot legend the first simulation
            p1 = stairs(neur_act.ISI,'LineWidth',1,'color',blueGRADIENTflexible(i,length(parameter_vec)));
            hold on;
        elseif i == length(parameter_vec)        % for inserting in the plot legend the last simulation
            p2 = stairs(neur_act.ISI,'LineWidth',1,'color',blueGRADIENTflexible(i,length(parameter_vec)));
            hold on;
        else
            stairs(neur_act.ISI,'LineWidth',1.5,'color',blueGRADIENTflexible(i,length(parameter_vec)));
            hold on;
        end

    end
    
    title('ISI of membrane potential');
    ylabel('ISI','FontName','Arial','FontSize',14);
    str1 = strcat(parameter_name,' = ', string(parameter_lower));
    str2 = strcat(parameter_name,' = ', string(parameter_upper));
    if filt
        str3 = strcat(model,' filtered');
    else
        str3 = model;
    end
    legend([p1(1), p2(1),p3(1)],{str1,str2,str3},'FontSize',14);                % shows the parameter range of the plot


    %% Plot of spike count
    flag = 1;
    
    % Sequence of simulations in the parameter range of interest
    for i=1:length(parameter_vec)

        if flag                 % Initialisation of figure for plot of membrane potential
            figID = figID+1;
            figure(figID);
            flag = 0;
        end

        [t,Y] = ode45(@(t,Y)Hindmarsh_Rose(t,Y,parameter_mat(:,i)),t_span,x_0,simoptions);
        [neur_act] = HR_neuronal_activity(Y(:,1),t,x_R);

        if i == 1                                % for inserting in the plot legend the first simulation
            p1 = stem(neur_act.spike_count,'LineWidth',1,'color',blueGRADIENTflexible(i,length(parameter_vec)));
            hold on;
        elseif i == length(parameter_vec)        % for inserting in the plot legend the last simulation
            p2 = stem(neur_act.spike_count,'LineWidth',1,'color',blueGRADIENTflexible(i,length(parameter_vec)));
            hold on;
        else
            stem(neur_act.spike_count,'LineWidth',1,'color',blueGRADIENTflexible(i,length(parameter_vec)));
            hold on;
        end

    end

    p3 = stem(neur_act_mean.spike_count,'r');
    p3.LineWidth = 2;
    
    title('Spike count of membrane potential');
    xlabel('Neuronal activity cycle','FontName','Arial','FontSize',14);
    ylabel('Spike count','FontName','Arial','FontSize',14);
    str1 = strcat(parameter_name,' = ', string(parameter_lower));
    str2 = strcat(parameter_name,' = ', string(parameter_upper));
    if filt
        str3 = strcat(model,' filtered');
    else
        str3 = model;
    end
    legend([p1(1), p2(1),p3(1)],{str1,str2,str3},'FontSize',14);                % shows the parameter range of the plot

    %% Plot of the duty cycle
    flag = 1;
    
    % Sequence of simulations in the parameter range of interest
    for i=1:length(parameter_vec)

        if flag                 % Initialisation of figure for plot of membrane potential
            figID = figID+1;
            figure(figID);
            flag = 0;
        end

        [t,Y] = ode45(@(t,Y)Hindmarsh_Rose(t,Y,parameter_mat(:,i)),t_span,x_0,simoptions);
        [neur_act] = HR_neuronal_activity(Y(:,1),t,x_R);

        if i == 1                                % for inserting in the plot legend the first simulation
            p1 = stem(100*neur_act.duty_cycle,'LineWidth',1,'color',blueGRADIENTflexible(i,length(parameter_vec)));
            hold on;
        elseif i == length(parameter_vec)        % for inserting in the plot legend the last simulation
            p2 = stem(100*neur_act.duty_cycle,'LineWidth',1,'color',blueGRADIENTflexible(i,length(parameter_vec)));
            hold on;
        else
            stem(100*neur_act.duty_cycle,'LineWidth',1,'color',blueGRADIENTflexible(i,length(parameter_vec)));
            hold on;
        end

    end

    p3 = stem(100*neur_act_mean.duty_cycle,'r');
    p3.LineWidth = 2;
    
    title('Duty cycle of membrane potential');
    xlabel('Neuronal activity cycle','FontName','Arial','FontSize',14);
    ylabel('Duty cycle [%]','FontName','Arial','FontSize',14);
    str1 = strcat(parameter_name,' = ', string(parameter_lower));
    str2 = strcat(parameter_name,' = ', string(parameter_upper));
    if filt
        str3 = strcat(model,' filtered');
    else
        str3 = model;
    end
    legend([p1(1), p2(1),p3(1)],{str1,str2,str3},'FontSize',14);                % shows the parameter range of the plot
    
    
end