function [results,tEnd] = surrogate_model(states,parameters,inputs,surr_model,simulation_opts,model)
% This function computes mean and variance based on the uncertainty
% structure specified by the structs "states" and "parameters" according to
% the surrogate model specified by "surr_model" and the simulation options
% specified by "simulation_opts". The time taken to compute the results is
% returned as well.

    % Setting of simulation parameters
    simoptions.tspan = [simulation_opts{1,1}, simulation_opts{1,2}];
    simoptions.dt = simulation_opts{1,3};
    simoptions.setup = odeset;
    simoptions.solver = simulation_opts{1,4};

    %% Run simulations

    if surr_model == 1  % Galerkin approach

        pce_order = simulation_opts{1,7};

        tStart = tic;

        % Compose the PCE system and write files
        sys = PoCETcompose(states,parameters,inputs,[],pce_order);                                 
        PoCETwriteFiles(sys,'deterministic_expanded_system',[],'nominal_system');

        % Run galerkin-PCE simulation
        results = PoCETsimGalerkin(sys,'deterministic_expanded_system',[],simoptions);  
        
        % Compute moments from simulation results
        sys.MomMats = PoCETmomentMatrices(sys,2);
        results.x.moments = PoCETcalcMoments(sys,sys.MomMats,results.x.pcvals);
        
        % Computation time
        tEnd = toc(tStart);

    elseif surr_model == 2  % Collocation approach
        
        pce_order = simulation_opts{1,7};
        colloc_samples = simulation_opts{1,6};

        tStart = tic;

        % Compose the PCE system and write files
        sys = PoCETcompose(states,parameters,inputs,[],pce_order);                                 
        if model ~= 2 && model ~= 4
            PoCETwriteFiles(sys,'deterministic_expanded_system',[],'nominal_system');
        end

        % Run collocation-PCE simulation
        basis = PoCETsample(sys,'basis',colloc_samples);    % 'basis' samples from the stochastic basis \xi and not the actual random variable
        if model == 2
            results = PoCETsimCollocation(sys,'Epileptor_nominal_system',[],basis,simoptions);  % problem with piecewise functions for f_1 and f_2
         elseif model == 4
             results = PoCETsimCollocation(sys,'autoactivating_feedback_nominal_system',[],basis,simoptions);   % I have no idea what the problem is with PoCETwriteFiles()
        else
            results = PoCETsimCollocation(sys,'nominal_system',[],basis,simoptions);
        end  
        
        % Compute moments from simulation results for quantities of interest
        sys.MomMats = PoCETmomentMatrices(sys,2);
        if model == 1
            results.x.moments = PoCETcalcMoments(sys,sys.MomMats,results.x.pcvals);
        elseif model == 2
            results.x_1.moments = PoCETcalcMoments(sys,sys.MomMats,results.x_1.pcvals);
            results.x_4.moments = PoCETcalcMoments(sys,sys.MomMats,results.x_4.pcvals);
        elseif model == 3
            results.y_2.moments = PoCETcalcMoments(sys,sys.MomMats,results.y_2.pcvals);
            results.y_3.moments = PoCETcalcMoments(sys,sys.MomMats,results.y_3.pcvals);
        elseif model == 4
            results.x.moments = PoCETcalcMoments(sys,sys.MomMats,results.x.pcvals);
        end

        % Computation time
        tEnd = toc(tStart);

    elseif surr_model == 3  % Monte Carlo

        mc_samples = simulation_opts{1,5};
        pce_order = simulation_opts{1,7};

        tStart = tic;

        % Compose the PCE system and write files
        sys = PoCETcompose(states,parameters,inputs,[],pce_order);
        if model ~= 2 && model ~= 4
            PoCETwriteFiles(sys,'deterministic_expanded_system',[],'nominal_system');
        end

        % Run Monte-Carlo simulations
        if model == 2
            results = PoCETsimMonteCarlo(sys,'Epileptor_nominal_system',[],mc_samples,simoptions,'method','moments');   % problem with piecewise functions for f_1 and f_2
        elseif model == 4
            results = PoCETsimMonteCarlo(sys,'autoactovating_feedback_nominal_system',[],mc_samples,simoptions,'method','moments'); % I have no idea what the problem is with PoCETwriteFiles()
        else
            results = PoCETsimMonteCarlo(sys,'nominal_system',[],mc_samples,simoptions,'method','moments');
        end

        % Computation time
        tEnd = toc(tStart);
    end

end