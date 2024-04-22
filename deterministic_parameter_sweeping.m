function figID = deterministic_parameter_sweeping(nom_parameters,x_0,uncertainty,simulation_opts, figID, model)
    % This function computes performances on an ensemble of simulations in
    % which one or more parameters vary in a range. In particular the
    % performances for the Hindmarsh-Rose model are the maximum spike 
    % amplitude evolution w.r.t. a parameter change, the average 
    % membrane potential evolution w.r.t. a parameter change, the 
    % variance of the membrane potential evolution w.r.t. a parameter 
    % change and evolution w.r.t. time of the membrane potential as one
    % parameter changes in a range.
    
    % Setting of simulation parameters
    parameter_number = uncertainty{1,1};
    interval = uncertainty{1,2};
    n_plots = simulation_opts{1,8};
    t_span = [simulation_opts{1,1}, simulation_opts{1,2}];
    simoptions.dt = simulation_opts{1,3};
    simoptions.setup = odeset;
    simoptions.solver = {1,4};

    % Name assignments for x-axis of the plot
    parameter_name = parameter_name_assignment(parameter_number,model);
   
    % Setting for plot of membrane potential evolution w.r.t. a changing parameter
    lightBLUE = [0.356862745098039,0.811764705882353,0.956862745098039];
    darkBLUE = [0.0196078431372549,0.0745098039215686,0.670588235294118];
    blueGRADIENTflexible = @(i,N) lightBLUE + (darkBLUE-lightBLUE)*((i-1)/(N-1));
    
    parameter = zeros(1,length(parameter_number));
    parameter_lower = parameter;
    parameter_upper = parameter;
    resolution = parameter;
    parameter_vec = cell(1,length(parameter_number));
    parameter_mat = cell(1,length(parameter_number));
    if model == 1
        result = cell(1,length(parameter_number));
        for k=1:length(parameter_number)
            result{1,k} = zeros(3,length(parameter_vec{1,k}));
        end
        
    end

    for k=1:length(parameter_number)
        % Creation of a vector of different values of a parameter of interest
        parameter(k) = nom_parameters(parameter_number(k)); 
        parameter_lower(k) = parameter(k)-interval(k)/2; 
        parameter_upper(k) = parameter(k)+interval(k)/2;
        resolution(k) = (parameter_upper(k)-parameter_lower(k))/n_plots;
        parameter_vec{1,k} = parameter_lower(k):resolution(k):parameter_upper(k);
        
        % Creation of a matrix of parameters in which only the parameter of interest changes
        parameter_mat{1,k} = kron(nom_parameters',ones(1,length(parameter_vec{1,k})));
        parameter_mat{1,k}(parameter_number(k),:) = parameter_vec{1,k};
    end

    flag = 1;
    
    for k=1:length(parameter_number)
        % Sequence of simulations in the parameter range of interest
        for i=1:length(parameter_vec{1,k})
    
            if flag                 % Initialisation of figure for plot of membrane potential
                figID = figID+1;
                figure(figID);
                flag = 0;
            end
            
            if model == 1
                [t,Y] = ode45(@(t,Y)Hindmarsh_Rose(t,Y,parameter_mat{1,k}(:,i)),t_span,x_0,simoptions);
            elseif model == 2
                [t,Y] = ode45(@(t,Y)Epileptor(t,Y,parameter_mat{1,k}(:,i)),t_span,x_0,simoptions);
            elseif model == 3
                [t,Y] = ode45(@(t,Y)Jansen_Rit(t,Y,parameter_mat{1,k}(:,i)),t_span,x_0,simoptions);
            elseif model == 4
                [t,Y] = ode45(@(t,Y)autoactivating_feedback_loop(t,Y,parameter_mat{1,k}(:,i)),t_span,x_0,simoptions);
            end
            
            if model == 1
                result{1,k}(1,i) = mean(maxk(Y(:,1),3));   % mean of the 3 maximum values is considered for robustness w.r.t. initial conditions
                result{1,k}(2,i) = mean(Y(:,1));
                result{1,k}(3,i) = var(Y(:,1));
            end
    
            if i == 1                                % for inserting in the plot legend the first simulation
                if model == 1
                    p1 = plot(t,Y(:,1),'LineWidth',1.5,'color',blueGRADIENTflexible(i,length(parameter_vec{1,k})));
                    hold on;
                    plot(t, nom_parameters(8)*ones(1,length(t)),'LineWidth',1.5);
                    hold on;
                elseif model == 2
                    p1 = plot(t,-Y(:,1)+Y(:,4),'LineWidth',1.5,'color',blueGRADIENTflexible(i,length(parameter_vec{1,k})));
                    hold on;
                elseif model == 3
                    p1 = plot(t,Y(:,2)-Y(:,3),'LineWidth',1.5,'color',blueGRADIENTflexible(i,length(parameter_vec{1,k})));
                    hold on;
                elseif model == 4
                    p1 = plot(t,Y,'LineWidth',1.5,'color',blueGRADIENTflexible(i,length(parameter_vec{1,k})));
                    hold on;
                end
            elseif i == length(parameter_vec{1,k})        % for inserting in the plot legend the last simulation
                if model == 1
                    p2 = plot(t,Y(:,1),'LineWidth',1.5,'color',blueGRADIENTflexible(i,length(parameter_vec{1,k})));
                    hold on;
                elseif model == 2
                    p2 = plot(t,-Y(:,1)+Y(:,4),'LineWidth',1.5,'color',blueGRADIENTflexible(i,length(parameter_vec{1,k})));
                    hold on;
                elseif model == 3
                    p2 = plot(t,Y(:,2)-Y(:,3),'LineWidth',1.5,'color',blueGRADIENTflexible(i,length(parameter_vec{1,k})));
                    hold on;
                elseif model == 4
                    p2 = plot(t,Y,'LineWidth',1.5,'color',blueGRADIENTflexible(i,length(parameter_vec{1,k})));
                    hold on;
                end
            else
                if model == 1
                    plot(t,Y(:,1),'LineWidth',1.5,'color',blueGRADIENTflexible(i,length(parameter_vec{1,k})));
                    hold on;
                elseif model == 2
                    plot(t,-Y(:,1)+Y(:,4),'LineWidth',1.5,'color',blueGRADIENTflexible(i,length(parameter_vec{1,k})));
                    hold on;
                elseif model == 3
                    plot(t,Y(:,2)-Y(:,3),'LineWidth',1.5,'color',blueGRADIENTflexible(i,length(parameter_vec{1,k})));
                    hold on;
                elseif model == 4
                    plot(t,Y,'LineWidth',1.5,'color',blueGRADIENTflexible(i,length(parameter_vec{1,k})));
                    hold on;
                end
            end
    
        end
        if model == 1
            title('Membrane potential');
            xlabel('time','FontName','Arial','FontSize',14);
            ylabel('x','FontName','Arial','FontSize',14);
        elseif model == 2
            xlabel('time','FontName','Arial','FontSize',14);
            ylabel('-x_1+x_4','FontName','Arial','FontSize',14);
        elseif model == 3
            xlabel('time','FontName','Arial','FontSize',14);
            ylabel('y_2-y_3','FontName','Arial','FontSize',14);
        elseif model == 4
            xlabel('time','FontName','Arial','FontSize',14);
            ylabel('x','FontName','Arial','FontSize',14);
        end
    
        str1 = strcat(parameter_name{1,k},' = ', string(parameter_lower(k)));
        str2 = strcat(parameter_name{1,k},' = ', string(parameter_upper(k)));
        legend([p1(1), p2(1)],{str1,str2},'FontSize',14);   % shows the parameter range of the plot
        
        if model == 1
            % Plot the spike amplitude evolution
            figID = figID+1;
            figure(figID);
            plot(parameter_vec{1,k},result{1,k}(1,:),'LineWidth',1.5);
            xlabel(parameter_name{1,k},'FontName','Arial','FontSize',14);  
            ylabel('Max spike amplitude','FontName','Arial','FontSize',14);
            str = strcat('Simulation time is',{' '},string(simulation_opts{1,2}));
            annotation('textbox',[0.2,0.5,0.3,0.3],'String',str,'FitBoxToText','on');
        
            % Plot the average membrane potential evolution
            figID = figID+1;
            figure(figID);
            plot(parameter_vec{1,k},result{1,k}(2,:),'LineWidth',1.5);
            % Plot of the line of the resting potential
            if parameter_number(k) ~= 8
                hold on;
                plot(parameter_vec{1,k}, nom_parameters(8)*ones(1,length(parameter_vec{1,k})),'LineWidth',1.5);
            end
            xlabel(parameter_name{1,k},'FontName','Arial','FontSize',14);  
            ylabel('Average spike amplitude','FontName','Arial','FontSize',14);
            str = strcat('Simulation time is',{' '},string(simulation_opts{1,2}));
            annotation('textbox',[0.2,0.5,0.3,0.3],'String',str,'FitBoxToText','on');
        
            % Plot the variance of the membrane potential evolution
            figID = figID+1;
            figure(figID);
            plot(parameter_vec{1,k},result{1,k}(1,:),'LineWidth',1.5);
            xlabel(parameter_name{1,k},'FontName','Arial','FontSize',14);  
            ylabel('Variance spike amplitude','FontName','Arial','FontSize',14);
            str = strcat('Simulation time is',{' '},string(simulation_opts{1,2}));
            annotation('textbox',[0.2,0.5,0.3,0.3],'String',str,'FitBoxToText','on');
        end
    
        flag = 1;

    end
end