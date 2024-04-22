function [figID,disagreement_times_vec,results1,results2] = UQ_comparison(states,parameters,inputs,ground_truth,comparison,simulation_opts,figID)
% This function compares the first two moments computed by two different
% approaches: either Galerkin approach vs Monte Carlo, Collocation approach
% vs Monte Carlo, or Galerkin approach vs Collocation approach. For
% comparison, a ground truth is considered with a confidence interval to
% show when the two approaches violate the confidence band.

    % Setting of simulation parameters
    simoptions.tspan = [simulation_opts{1,1}, simulation_opts{1,2}];
    simoptions.dt = simulation_opts{1,3};
    simoptions.setup = odeset;
    simoptions.solver = simulation_opts{1,4};
    threshold_abs = simulation_opts{1,9};
    threshold_rel = simulation_opts{1,10};

    disagreement_times_vec = zeros(1,8);

    % Resize of the ground truth according to the final simulation time
    ground_truth.time = ground_truth.time(:,ground_truth.time <= simulation_opts{1,2});

    %% Run simulations

    if comparison == 1  % Galerkin approach vs Monte Carlo

        pce_order = simulation_opts{1,7};
        mc_samples = simulation_opts{1,5};

        tStart1 = tic;

        % Compose the PCE system and write files for galerkin-PCE and MC simulation
        sys = PoCETcompose(states,parameters,inputs,[],pce_order);                                 
        PoCETwriteFiles(sys,'deterministic_expanded_system',[],'nominal_system');
        % Run galerkin-PCE simulation
        results1 = PoCETsimGalerkin(sys,'deterministic_expanded_system',[],simoptions);  
        
        % Compute moments from simulation results
        sys.MomMats = PoCETmomentMatrices(sys,2);
        results1.x.moments = PoCETcalcMoments(sys,sys.MomMats,results1.x.pcvals);

        tEnd1 = toc(tStart1);

        tStart2 = tic;
            
        % Run Monte-Carlo simulations
        results2 = PoCETsimMonteCarlo(sys,'nominal_system',[],mc_samples,simoptions,'method','moments');

        tEnd2 = toc(tStart2);

    elseif comparison == 2  % Collocation approach vs Monte Carlo

        pce_order = simulation_opts{1,7};
        colloc_samples = simulation_opts{1,6};
        mc_samples = simulation_opts{1,5};

        tStart1 = tic;

        % Compose the PCE system and write files for collocation-PCE and MC simulation
        sys = PoCETcompose(states,parameters,inputs,[],pce_order);                                 
        PoCETwriteFiles(sys,'deterministic_expanded_system',[],'nominal_system');
        % Run collocation-PCE simulation
        basis = PoCETsample(sys,'basis',colloc_samples);    % 'basis' samples from the stochastic basis \xi and not the actual random variable
        results1 = PoCETsimCollocation(sys,'nominal_system',[],basis,simoptions);  
        
        % Compute moments from simulation results
        sys.MomMats = PoCETmomentMatrices(sys,2);
        results1.x.moments = PoCETcalcMoments(sys,sys.MomMats,results1.x.pcvals);

        tEnd1 = toc(tStart1);

        tStart2 = tic;
            
        % Run Monte-Carlo simulations
        results2 = PoCETsimMonteCarlo(sys,'nominal_system',[],mc_samples,simoptions,'method','moments');

        tEnd2 = toc(tStart2);

    elseif comparison == 3  % Galerkin approach vs collocation approach

        pce_order = simulation_opts{1,7};
        colloc_samples = simulation_opts{1,6};

        tStart1 = tic;

        % Compose the PCE system and write files for galerkin-PCE
        sys = PoCETcompose(states,parameters,inputs,[],pce_order);                                 
        PoCETwriteFiles(sys,'deterministic_expanded_system',[],'nominal_system');
        % Run galerkin-PCE simulation
        results1 = PoCETsimGalerkin(sys,'deterministic_expanded_system',[],simoptions);  
        
        % Compute moments from simulation results
        sys.MomMats = PoCETmomentMatrices(sys,2);
        results1.x.moments = PoCETcalcMoments(sys,sys.MomMats,results1.x.pcvals);

        tEnd1 = toc(tStart1);

        tStart2 = tic;

        % Run collocation-PCE simulation
        basis = PoCETsample(sys,'basis',colloc_samples);    % 'basis' samples from the stochastic basis \xi and not the actual random variable
        results2 = PoCETsimCollocation(sys,'nominal_system',[],basis,simoptions);  
        
        % Compute moments from simulation results
        results2.x.moments = PoCETcalcMoments(sys,sys.MomMats,results2.x.pcvals);

        tEnd2 = toc(tStart2);

    end

    %% Compute the disagreement times
    abs = 1;    % option for absolute error
    rel = 2;    % option for relative error

    [time1_mean_abs,time1_var_abs,time2_mean_abs,time2_var_abs] = disagreement_times(ground_truth,results1,results2,abs,threshold_abs);
    disagreement_times_vec(1) = time1_mean_abs;
    disagreement_times_vec(2) = time1_var_abs;
    disagreement_times_vec(3) = time2_mean_abs;
    disagreement_times_vec(4) = time2_var_abs;
    
    [time1_mean_rel,time1_var_rel,time2_mean_rel,time2_var_rel] = disagreement_times(ground_truth,results1,results2,rel,threshold_rel);
    disagreement_times_vec(5) = time1_mean_rel;
    disagreement_times_vec(6) = time1_var_rel;
    disagreement_times_vec(7) = time2_mean_rel;
    disagreement_times_vec(8) = time2_var_rel;


    %% Plot results for the mean with absolute error
    figID = figID+1;
    figure(figID);

    p1 = plot(ground_truth.time,ground_truth.x.moments(1,1:length(ground_truth.time)),'r','LineWidth',2.5);
    hold on;
    % Confidence bounds to detect disagreement
    plot(ground_truth.time,ground_truth.x.moments(1,1:length(ground_truth.time))+threshold_abs*ones(1,length(ground_truth.time)),'--r','LineWidth',2.5);    
    hold on;
    plot(ground_truth.time,ground_truth.x.moments(1,1:length(ground_truth.time))-threshold_abs*ones(1,length(ground_truth.time)),'--r','LineWidth',2.5);
    hold on;

    p2 = plot(results1.time,results1.x.moments(1,:),'g','LineWidth',1.5);
    hold on;

    p3 = plot(results2.time,results2.x.moments(1,:),'b','LineWidth',1.5);
    hold on;

    xline(time1_mean_abs,'g:','LineWidth',1.5);
    hold on;
    xline(time2_mean_abs,'b:','LineWidth',1.5);

    title('Mean with absolute admissible error');
    xlabel('time','FontName','Arial','FontSize',14);
    ylabel('Mean x','FontName','Arial','FontSize',14);

    if comparison == 1      % Galerkin approach vs Monte Carlo
        legend([p1(1), p2(1),p3(1)],{'Ground truth','Galerkin gPC','Monte Carlo'},'FontSize',14);
    elseif comparison == 2  % Collocation approach vs Monte Carlo
        legend([p1(1), p2(1),p3(1)],{'Ground truth','Collocation gPC','Monte Carlo'},'FontSize',14);
    elseif comparison == 3  % Galerkin approach vs collocation approach
        legend([p1(1), p2(1),p3(1)],{'Ground truth','Galerkin gPC','Collocation gPC'},'FontSize',14);
    end

    %% Plot results for the mean with relative error
    figID = figID+1;
    figure(figID);

    p1 = plot(ground_truth.time,ground_truth.x.moments(1,1:length(ground_truth.time)),'r','LineWidth',2.5);
    hold on;
    % Confidence bounds to detect disagreement
    plot(ground_truth.time,(1+threshold_rel)*ground_truth.x.moments(1,1:length(ground_truth.time)),'--r','LineWidth',2.5);    
    hold on;
    plot(ground_truth.time,(1-threshold_rel)*ground_truth.x.moments(1,1:length(ground_truth.time)),'--r','LineWidth',2.5);
    hold on;

    p2 = plot(results1.time,results1.x.moments(1,:),'g','LineWidth',1.5);
    hold on;

    p3 = plot(results2.time,results2.x.moments(1,:),'b','LineWidth',1.5);
    hold on;

    xline(time1_mean_rel,'g:','LineWidth',1.5);
    hold on;
    xline(time2_mean_rel,'b:','LineWidth',1.5);

    title('Mean with relative admissible error');
    xlabel('time','FontName','Arial','FontSize',14);
    ylabel('Mean x','FontName','Arial','FontSize',14);

    if comparison == 1      % Galerkin approach vs Monte Carlo
        legend([p1(1), p2(1),p3(1)],{'Ground truth','Galerkin gPC','Monte Carlo'},'FontSize',14);
    elseif comparison == 2  % Collocation approach vs Monte Carlo
        legend([p1(1), p2(1),p3(1)],{'Ground truth','Collocation gPC','Monte Carlo'},'FontSize',14);
    elseif comparison == 3  % Galerkin approach vs collocation approach
        legend([p1(1), p2(1),p3(1)],{'Ground truth','Galerkin gPC','Collocation gPC'},'FontSize',14);
    end


    %% Plot results for the variance with absolute error
    figID = figID+1;
    figure(figID);

    p1 = plot(ground_truth.time,ground_truth.x.moments(2,1:length(ground_truth.time)),'r','LineWidth',2.5);
    hold on;
    % Confidence bounds to detect disagreement
    plot(ground_truth.time,ground_truth.x.moments(2,1:length(ground_truth.time))+threshold_abs*ones(1,length(ground_truth.time)),'--r','LineWidth',2.5);    
    hold on;
    plot(ground_truth.time,ground_truth.x.moments(2,1:length(ground_truth.time))-threshold_abs*ones(1,length(ground_truth.time)),'--r','LineWidth',2.5);
    hold on;

    p2 = plot(results1.time,results1.x.moments(2,:),'g','LineWidth',1.5);
    hold on;

    p3 = plot(results2.time,results2.x.moments(2,:),'b','LineWidth',1.5);
    hold on;

    xline(time1_var_abs,'g:','LineWidth',1.5);
    hold on;
    xline(time2_var_abs,'b:','LineWidth',1.5);

    title('Variance with absolute admissible error');
    xlabel('time','FontName','Arial','FontSize',14);
    ylabel('Variance x','FontName','Arial','FontSize',14);

    if comparison == 1      % Galerkin approach vs Monte Carlo
        legend([p1(1), p2(1),p3(1)],{'Ground truth','Galerkin gPC','Monte Carlo'},'FontSize',14);
    elseif comparison == 2  % Collocation approach vs Monte Carlo
        legend([p1(1), p2(1),p3(1)],{'Ground truth','Collocation gPC','Monte Carlo'},'FontSize',14);
    elseif comparison == 3  % Galerkin approach vs collocation approach
        legend([p1(1), p2(1),p3(1)],{'Ground truth','Galerkin gPC','Collocation gPC'},'FontSize',14);
    end

    %% Plot results for the variance with relative error
    figID = figID+1;
    figure(figID);

    p1 = plot(ground_truth.time,ground_truth.x.moments(2,1:length(ground_truth.time)),'r','LineWidth',2.5);
    hold on;
    % Confidence bounds to detect disagreement
    plot(ground_truth.time,(1+threshold_rel)*ground_truth.x.moments(2,1:length(ground_truth.time)),'--r','LineWidth',2.5);    
    hold on;
    plot(ground_truth.time,(1-threshold_rel)*ground_truth.x.moments(2,1:length(ground_truth.time)),'--r','LineWidth',2.5);
    hold on;

    p2 = plot(results1.time,results1.x.moments(2,:),'g','LineWidth',1.5);
    hold on;

    p3 = plot(results2.time,results2.x.moments(2,:),'b','LineWidth',1.5);
    hold on;

    xline(time1_var_rel,'g:','LineWidth',1.5);
    hold on;
    xline(time2_var_rel,'b:','LineWidth',1.5);
    
    title('Variance with relative admissible error');
    xlabel('time','FontName','Arial','FontSize',14);
    ylabel('Variance x','FontName','Arial','FontSize',14);

    if comparison == 1      % Galerkin approach vs Monte Carlo
        legend([p1(1), p2(1),p3(1)],{'Ground truth','Galerkin gPC','Monte Carlo'},'FontSize',14);
    elseif comparison == 2  % Collocation approach vs Monte Carlo
        legend([p1(1), p2(1),p3(1)],{'Ground truth','Collocation gPC','Monte Carlo'},'FontSize',14);
    elseif comparison == 3  % Galerkin approach vs collocation approach
        legend([p1(1), p2(1),p3(1)],{'Ground truth','Galerkin gPC','Collocation gPC'},'FontSize',14);
    end
    
    % Output execution times
    if comparison == 1      % Galerkin approach vs Monte Carlo
        fprintf('\n');
        disp(['Galerkin approach expansion order ',num2str(pce_order),' took ',num2str(tEnd1), ' seconds to execute.']);
        disp(['Monte Carlo approach with ',num2str(mc_samples),' samples took ',num2str(tEnd2), ' seconds to execute.']);
        fprintf('\n');
    elseif comparison == 2  % Collocation approach vs Monte Carlo
        fprintf('\n');
        disp(['Collocation approach expansion order ',num2str(pce_order),' with ',num2str(colloc_samples),' samples took ',num2str(tEnd1), ' seconds to execute.']);
        disp(['Monte Carlo approach with ',num2str(mc_samples),' samples took ',num2str(tEnd2), ' seconds to execute.']);
        fprintf('\n');
    elseif comparison == 3  % Galerkin approach vs collocation approach
        fprintf('\n');
        disp(['Galerkin approach expansion order ',num2str(pce_order),' took ',num2str(tEnd1), ' seconds to execute.']);
        disp(['Collocation approach expansion order ',num2str(pce_order),' with ',num2str(colloc_samples),' samples took ',num2str(tEnd2), ' seconds to execute.']);
        fprintf('\n');
    end

end
