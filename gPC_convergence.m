function figID = gPC_convergence(states,parameters,inputs,ground_truth,max_order,intrusive,simulation_opts,figID)
    % This function computes the first two moments of the gPC expansion for
    % increasing expansion order from 1 to max_order. The figures contain
    % mean and variance of all the expansion orders, the absolute value
    % of the difference of the two moments w.r.t. the ground truth and the
    % relative error.
    
    % Setting of simulation parameters and plot parameters
    simoptions.tspan = [simulation_opts{1,1}, simulation_opts{1,2}];
    simoptions.dt = simulation_opts{1,3};
    simoptions.setup = odeset;
    simoptions.solver = simulation_opts{1,4};

    if intrusive == 2
        colloc_samples = simulation_opts{1,6};
    end

    % Resize of the ground truth according to the final simulation time
    ground_truth.time = ground_truth.time(:,ground_truth.time <= simulation_opts{1,2});
    
    % Settings for figures and plot of the ground truth
    fig_mean = figID+1;
    fig_mean_abs_err = figID+2;
    fig_mean_rel_err = figID+3;
    fig_variance = figID+4;
    fig_variance_abs_err = figID+5;
    fig_variance_rel_err = figID+6;

    hplot1 = cell(1,6);
    hplot2 = cell(1,5);
    hplot3 = cell(1,5);
    hplot4 = cell(1,6);
    hplot5 = cell(1,5);
    hplot6 = cell(1,5);
    
    figure(fig_mean);
    hplot1{1,6} = plot(ground_truth.time,ground_truth.x.moments(1,1:length(ground_truth.time)),'r','LineWidth',2.5);
    hold on;
    title('Mean');
    xlabel('time','FontName','Arial','FontSize',14);
    ylabel('Mean x','FontName','Arial','FontSize',14);

    set(gca,'ColorOrderIndex',1);

    figure(fig_variance);
    hplot4{1,6} = plot(ground_truth.time,ground_truth.x.moments(2,1:length(ground_truth.time)),'r','LineWidth',2.5);
    hold on;
    title('Variance');
    xlabel('time','FontName','Arial','FontSize',14);
    ylabel('Variance x','FontName','Arial','FontSize',14);

    set(gca,'ColorOrderIndex',1);
    
    % Expansion order loop
    for pce_order=1:max_order

        if intrusive == 1 % Galerkin approach

            tStart = tic;

            % Compose the PCE system and write files for galerkin-PCE
            sys = PoCETcompose(states,parameters,inputs,[],pce_order);                                 
            PoCETwriteFiles(sys,'deterministic_expanded_system',[],'nominal_system');
            % Run galerkin-PCE simulation
            results = PoCETsimGalerkin(sys,'deterministic_expanded_system',[],simoptions);   % contains solution of the expanded variables
                
            % Compute moments from simulation results
            sys.MomMats = PoCETmomentMatrices(sys,2);
            results.x.moments = PoCETcalcMoments(sys,sys.MomMats,results.x.pcvals);

            tEnd = toc(tStart);

        elseif intrusive == 2 % Collocation approach

            tStart = tic;

            % Compose the PCE system and write files for collocation-PCE
            sys = PoCETcompose(states,parameters,inputs,[],pce_order);                                 
            PoCETwriteFiles(sys,'deterministic_expanded_system',[],'nominal_system');
            % Run collocation-PCE simulation
            basis = PoCETsample(sys,'basis',colloc_samples);    % 'basis' samples from the stochastic basis \xi and not the actual random variable
            results = PoCETsimCollocation(sys,'nominal_system',[],basis,simoptions);

            % Compute moments from simulation results
            sys.MomMats = PoCETmomentMatrices(sys,2);
            results.x.moments = PoCETcalcMoments(sys,sys.MomMats,results.x.pcvals);

            tEnd = toc(tStart);

        end

        % Plot mean of simulation
        figure(fig_mean);
        hplot1{1,mod(pce_order,5)+1} = plot(results.time,results.x.moments(1,:),'LineWidth',1.5);
        hold on;

        % Plot absolute error
        figure(fig_mean_abs_err);
        hplot2{1,mod(pce_order,5)+1} = plot(results.time,abs(results.x.moments(1,:)-ground_truth.x.moments(1,1:length(ground_truth.time))),'LineWidth',1.5);
        hold on;
        title('Absolute error mean');
        xlabel('time','FontName','Arial','FontSize',14);
        ylabel('Absolute error mean x','FontName','Arial','FontSize',14);
        
        % Plot relative error
        figure(fig_mean_rel_err);
        hplot3{1,mod(pce_order,5)+1} = plot(results.time,abs(results.x.moments(1,:)-ground_truth.x.moments(1,1:length(ground_truth.time)))./abs(ground_truth.x.moments(1,1:length(ground_truth.time))),'LineWidth',1.5);
        hold on;
        title('Relative error mean');
        xlabel('time','FontName','Arial','FontSize',14);
        ylabel('Relative error mean x','FontName','Arial','FontSize',14);

        % Plot variance of simulation
        figure(fig_variance);
        hplot4{1,mod(pce_order,5)+1} = plot(results.time,results.x.moments(2,:),'LineWidth',1.5);
        hold on;

        % Plot absolute error 
        figure(fig_variance_abs_err);
        hplot5{1,mod(pce_order,5)+1} = plot(results.time,abs(results.x.moments(2,:)-ground_truth.x.moments(2,1:length(ground_truth.time))),'LineWidth',1.5);
        hold on;
        title('Absolute error variance');
        xlabel('time','FontName','Arial','FontSize',14);
        ylabel('Absolute error variance x','FontName','Arial','FontSize',14);

        % Plot relative error 
        figure(fig_variance_rel_err);
        hplot6{1,mod(pce_order,5)+1} = plot(results.time,abs(results.x.moments(2,:)-ground_truth.x.moments(2,1:length(ground_truth.time)))./abs(ground_truth.x.moments(2,1:length(ground_truth.time))),'LineWidth',1.5);
        hold on;
        title('Relative error variance');
        xlabel('time','FontName','Arial','FontSize',14);
        ylabel('Relative error variance x','FontName','Arial','FontSize',14);

        % Program stops at each iteration, so that the user can decide when
        % to continue. Also execution time is printed.
        if pce_order < max_order
            if intrusive == 1
                fprintf('\n');
                disp(['Galerkin approach expansion order ',num2str(pce_order),' took ',num2str(tEnd), ' seconds to execute.']);
                fprintf('\n');
            elseif intrusive == 2
                fprintf('\n');
                disp(['Collocation approach expansion order ',num2str(pce_order),' with ',num2str(colloc_samples),' samples took ',num2str(tEnd), ' seconds to execute.']);
                fprintf('\n');
            end
            disp('Press any key to continue.');
            pause;
            % Keep only 5 graphs at a time for clarity, at each iteration
            % the oldest graph is deleted
            if mod(pce_order,5)+2 == 6
                delete(hplot1{1,1});
                delete(hplot2{1,1});
                delete(hplot3{1,1});
                delete(hplot4{1,1});
                delete(hplot5{1,1});
                delete(hplot6{1,1});
            else
                delete(hplot1{1,mod(pce_order,5)+2});
                delete(hplot2{1,mod(pce_order,5)+2});
                delete(hplot3{1,mod(pce_order,5)+2});
                delete(hplot4{1,mod(pce_order,5)+2});
                delete(hplot5{1,mod(pce_order,5)+2});
                delete(hplot6{1,mod(pce_order,5)+2});
            end
        else
            disp('Max expansion order reached.');
        end

    end

    figID = figID+6;
end