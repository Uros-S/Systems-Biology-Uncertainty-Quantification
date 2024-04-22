function [states,parameters,inputs] = HS_struct_create(nom_parameters,x_0,uncertainty)
    % This function creates the structs needed for the PoCET toolbox for
    % gPC expansions. First the nominal system is created and then uniform
    % uncertainty is added having parameter_of_interest at the center of an
    % interval with length uncertainty_delta. 

    parameters_of_interest = uncertainty{1,1};
    uncertainty_delta = uncertainty{1,2};
    
    % Creation of structs for nominal system
    states(1).name = 'x_1';
    states(1).dist = 'none';
    states(1).data = x_0(1);
    states(1).rhs = 'k_p*x_2*x_12 - k_pp*x_1*exp(-9.2369+1.591*Ea)*exp(-Ea/(8.3144598*(273.15+T_C)))';
            
    states(2).name = 'x_2';
    states(2).dist = 'none';
    states(2).data = x_0(2);
    states(2).rhs = '-k_p*x_2*x_12 + k_pp*x_1*exp(-9.2369+1.591*Ea)*exp(-Ea/(8.3144598*(273.15+T_C)))';
            
    states(3).name = 'x_3';
    states(3).dist = 'none';
    states(3).data = x_0(3);
    states(3).rhs = 'k_s*x_4 - k_sp*x_3*(x_2^5/(600+x_2^5))';

    states(4).name = 'x_4';
    states(4).dist = 'none';
    states(4).data = x_0(4);
    states(4).rhs = '-k_s*x_4 + k_sp*x_3*(x_2^5/(600+x_2^5))';
            
    states(5).name = 'x_5';
    states(5).dist = 'none';
    states(5).data = x_0(5);
    states(5).rhs = 'k_F*x_6 + k_piF*x_10 + k_FGp*x_9 - k_FG*x_7*x_5 - k_Fp*x_5*(x_4/(1+x_4)) - d_F*x_5';
            
    states(6).name = 'x_6';
    states(6).dist = 'none';
    states(6).data = x_0(6);
    states(6).rhs = '-k_F*x_6 + k_Fp*x_5*(x_4/(1+x_4)) + k_FsGp*x_8 - k_FsG*x_7*x_6';

    states(7).name = 'x_7';
    states(7).dist = 'none';
    states(7).data = x_0(7);
    states(7).rhs = 'k_FsGp*x_8 + k_FGp*x_9 - k_FG*x_5*x_7 - k_FsG*x_6*x_7';
            
    states(8).name = 'x_8';
    states(8).dist = 'none';
    states(8).data = x_0(8);
    states(8).rhs = 'k_FsG*x_6*x_7 + k_Fs*x_9 - k_Fsp*x_8 - k_FsGp*x_8';
            
    states(9).name = 'x_9';
    states(9).dist = 'none';
    states(9).data = x_0(9);
    states(9).rhs = 'k_Fsp*x_8 + k_FG*x_5*x_7 - k_FGp*x_9 - k_Fs*x_9';

    states(10).name = 'x_10';
    states(10).dist = 'none';
    states(10).data = x_0(10);
    states(10).rhs = 'k_piRF*x_8 + 0.02125*(d_RP*d_HP)/(k_piHP) - d_RF*x_10';
            
    states(11).name = 'x_11';
    states(11).dist = 'none';
    states(11).data = x_0(11);
    states(11).rhs = 'k_piRH*x_8 + 17.5*(d_RF*d_F)/(k_piF) - d_RHP*x_11';
            
    states(12).name = 'x_12';
    states(12).dist = 'none';
    states(12).data = x_0(12);
    states(12).rhs = 'k_piHP*x_11 - d_HP*x_12';


            
    parameters(1).name = 'k_p';
    parameters(1).dist = 'none';
    parameters(1).data = nom_parameters(1);
            
    parameters(2).name = 'k_pp';
    parameters(2).dist = 'none';
    parameters(2).data = nom_parameters(2);
            
    parameters(3).name = 'k_s';
    parameters(3).dist = 'none';
    parameters(3).data = nom_parameters(3);
            
    parameters(4).name = 'k_sp';
    parameters(4).dist = 'none';
    parameters(4).data = nom_parameters(4);
            
    parameters(5).name = 'k_F';
    parameters(5).dist = 'none';
    parameters(5).data = nom_parameters(5);
            
    parameters(6).name = 'k_Fp';
    parameters(6).dist = 'none';
    parameters(6).data = nom_parameters(6);
            
    parameters(7).name = 'k_FG';
    parameters(7).dist = 'none';
    parameters(7).data = nom_parameters(7);
            
    parameters(8).name = 'k_FGp';
    parameters(8).dist = 'none';
    parameters(8).data = nom_parameters(8);

    parameters(9).name = 'k_FsG';
    parameters(9).dist = 'none';
    parameters(9).data = nom_parameters(9);
            
    parameters(10).name = 'k_FsGp';
    parameters(10).dist = 'none';
    parameters(10).data = nom_parameters(10);
            
    parameters(11).name = 'k_Fs';
    parameters(11).dist = 'none';
    parameters(11).data = nom_parameters(11);
            
    parameters(12).name = 'k_Fsp';
    parameters(12).dist = 'none';
    parameters(12).data = nom_parameters(12);
            
    parameters(13).name = 'k_piRF';
    parameters(13).dist = 'none';
    parameters(13).data = nom_parameters(13);
            
    parameters(14).name = 'k_piRH';
    parameters(14).dist = 'none';
    parameters(14).data = nom_parameters(14);
            
    parameters(15).name = 'k_piHP';
    parameters(15).dist = 'none';
    parameters(15).data = nom_parameters(15);
            
    parameters(16).name = 'k_piF';
    parameters(16).dist = 'none';
    parameters(16).data = nom_parameters(16);

    parameters(17).name = 'd_F';
    parameters(17).dist = 'none';
    parameters(17).data = nom_parameters(17);
            
    parameters(18).name = 'd_HP';
    parameters(18).dist = 'none';
    parameters(18).data = nom_parameters(18);
            
    parameters(19).name = 'd_RF';
    parameters(19).dist = 'none';
    parameters(19).data = nom_parameters(19);
            
    parameters(20).name = 'd_RP';
    parameters(20).dist = 'none';
    parameters(20).data = nom_parameters(20);

    parameters(21).name = 'Ea';
    parameters(21).dist = 'none';
    parameters(21).data = nom_parameters(21);

    parameters(22).name = 'T_C';
    parameters(22).dist = 'none';
    parameters(22).data = nom_parameters(22);
    
    
    inputs = [];

    % Addition of uncertainty on parameter of interest
    for i=1:length(parameters_of_interest)
        parameters(parameters_of_interest(i)).dist = 'uniform';
        parameters(parameters_of_interest(i)).data = [nom_parameters(parameters_of_interest(i))-uncertainty_delta(i)/2, ...
                                                  nom_parameters(parameters_of_interest(i))+uncertainty_delta(i)/2];
    end

end