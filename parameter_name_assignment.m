function parameter_name = parameter_name_assignment(parameter_of_interest, model)
% This function outputs the string of the name of the parameter in position 
% parameter_number according to the model given by the number model.
    parameter_name = cell(1,length(parameter_of_interest));

    for i=1:length(parameter_of_interest)
        if model == 1
            if parameter_of_interest(i) == 1
                parameter_name{1,i} = 'a';
            elseif parameter_of_interest(i) == 2
                parameter_name{1,i} = 'b';
            elseif parameter_of_interest(i) == 3
                parameter_name{1,i} = 'c';
            elseif parameter_of_interest(i) == 4
                parameter_name{1,i} = 'd';
            elseif parameter_of_interest(i) == 5
                parameter_name{1,i} = 'I';
            elseif parameter_of_interest(i) == 6
                parameter_name{1,i} = 'r';
            elseif parameter_of_interest(i) == 7
                parameter_name{1,i} = 's';
            else
                parameter_name{1,i} = 'x_R';
            end
        elseif model == 2
            if parameter_of_interest(i) == 1
                parameter_name{1,i} = 'r_1';
            elseif parameter_of_interest(i) == 2
                parameter_name{1,i} = 'r_2';
            elseif parameter_of_interest(i) == 3
                parameter_name{1,i} = 'I_1';
            elseif parameter_of_interest(i) == 4
                parameter_name{1,i} = 'I_2';
            elseif parameter_of_interest(i) == 5
                parameter_name{1,i} = '\tau_0';
            elseif parameter_of_interest(i) == 6
                parameter_name{1,i} = 'm';
            elseif parameter_of_interest(i) == 7
                parameter_name{1,i} = '\tau_2';
            else
                parameter_name{1,i} = '\gamma';
            end
        elseif model == 3
            if parameter_of_interest(i) == 1
                parameter_name{1,i} = 'C';
            elseif parameter_of_interest(i) == 2
                parameter_name{1,i} = 'A';
            elseif parameter_of_interest(i) == 3
                parameter_name{1,i} = 'B';
            elseif parameter_of_interest(i) == 4
                parameter_name{1,i} = 'a';
            elseif parameter_of_interest(i) == 5
                parameter_name{1,i} = 'b';
            elseif parameter_of_interest(i) == 6
                parameter_name{1,i} = 'V_0';
            elseif parameter_of_interest(i) == 7
                parameter_name{1,i} = '\nu_{max}';
            elseif  parameter_of_interest(i) == 8
                parameter_name{1,i} = 'r';
            else
                parameter_name{1,i} = 'p';
            end
        elseif model == 4
            if parameter_of_interest(i) == 1
                parameter_name{1,i} = 'n';
            elseif parameter_of_interest(i) == 2
                parameter_name{1,i} = 'K';
            else
                parameter_name{1,i} = 'a';
            end
        elseif model == 5
            if parameter_of_interest(i) == 1
                parameter_name{1,i} = 'k_p';
            elseif parameter_of_interest(i) == 2
                parameter_name{1,i} = 'k_p''';
            elseif parameter_of_interest(i) == 3
                parameter_name{1,i} = 'k_s';
            elseif parameter_of_interest(i) == 4
                parameter_name{1,i} = 'k_s''';
            elseif parameter_of_interest(i) == 5
                parameter_name{1,i} = 'k_F';
            elseif parameter_of_interest(i) == 6
                parameter_name{1,i} = 'k_F''';
            elseif parameter_of_interest(i) == 7
                parameter_name{1,i} = 'k_{FG}';
            elseif parameter_of_interest(i) == 8
                parameter_name{1,i} = 'k_{FG}''';
            elseif parameter_of_interest(i) == 9
                parameter_name{1,i} = 'k_{F*G}';
            elseif parameter_of_interest(i) == 10
                parameter_name{1,i} = 'k_{F*G}''';
            elseif parameter_of_interest(i) == 11
                parameter_name{1,i} = 'k_{F*}';
            elseif parameter_of_interest(i) == 12
                parameter_name{1,i} = 'k_{F*}''';
            elseif parameter_of_interest(i) == 13
                parameter_name{1,i} = 'k_{\pi RF}';
            elseif parameter_of_interest(i) == 14
                parameter_name{1,i} = 'k_{\pi RH}';
            elseif parameter_of_interest(i) == 15
                parameter_name{1,i} = 'k_{\pi HP}';
            elseif parameter_of_interest(i) == 16
                parameter_name{1,i} = 'k_{\pi F}';
            elseif parameter_of_interest(i) == 17
                parameter_name{1,i} = 'd_F''';
            elseif parameter_of_interest(i) == 18
                parameter_name{1,i} = 'd_{HP}';
            elseif  parameter_of_interest(i) == 19
                parameter_name{1,i} = 'd_{RF}';
            elseif parameter_of_interest(i) == 20
                parameter_name{1,i} = 'd_{RP}';
            elseif parameter_of_interest(i) == 21
                parameter_name{1,i} = 'Ea';
            else
                parameter_name{1,i} = 'T_C';
            end
        end

    end
end