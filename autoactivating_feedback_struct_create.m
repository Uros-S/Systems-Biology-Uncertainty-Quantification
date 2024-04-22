function [states,parameters,inputs] = autoactivating_feedback_struct_create(nom_parameters,x_0,uncertainty)
    % This function creates the structs needed for the PoCET toolbox for
    % gPC expansions. First the nominal system is created and then uniform
    % uncertainty is added having parameter_of_interest at the center of an
    % interval with length uncertainty_delta. 

    parameters_of_interest = uncertainty{1,1};
    uncertainty_delta = uncertainty{1,2};
    
    % Creation of structs for nominal system
    states(1).name = 'x';
    states(1).dist = 'none';
    states(1).data = x_0;
    states(1).rhs = 'K + a*(x^n)/(1+x^n) - x';


            
    parameters(1).name = 'n';
    parameters(1).dist = 'none';
    parameters(1).data = nom_parameters(1);
            
    parameters(2).name = 'K';
    parameters(2).dist = 'none';
    parameters(2).data = nom_parameters(2);
            
    parameters(3).name = 'a';
    parameters(3).dist = 'none';
    parameters(3).data = nom_parameters(3);
    
    

    inputs = [];

    % Addition of uncertainty on parameter of interest
    for i=1:length(parameters_of_interest)
        parameters(parameters_of_interest(i)).dist = 'uniform';
        parameters(parameters_of_interest(i)).data = [nom_parameters(parameters_of_interest(i))-uncertainty_delta(i)/2, ...
                                                  nom_parameters(parameters_of_interest(i))+uncertainty_delta(i)/2];
    end

end