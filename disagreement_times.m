function [time1_mean,time1_var,time2_mean,time2_var] = disagreement_times(ground_truth,results1,results2,type_error,threshold)
% This function computes the disagreement times of results1 and results2
% w.r.t ground_truth considering both the mean and the variance. The type
% of confidence bound is specified by type_error and it is either an
% absolute error or a relative error. According to type_error a certain
% threshold is specified. If results do not violate the confidence bound
% then a NaN is returned.

    %% Find the disagreement times for results1 for the mean
    if type_error == 1
        temp_1 = (ground_truth.x.moments(1,1:length(ground_truth.time))+threshold*ones(1,length(ground_truth.time)))-results1.x.moments(1,:);
    elseif type_error == 2
        temp_1 = (1+threshold)*ground_truth.x.moments(1,1:length(ground_truth.time))-results1.x.moments(1,:);
    end
    
    indexOfMinDiff = find(temp_1 < 0, 1);
    time_1_a = results1.time(indexOfMinDiff);

    if type_error == 1
        temp_2 = results1.x.moments(1,:)-(ground_truth.x.moments(1,1:length(ground_truth.time))-threshold*ones(1,length(ground_truth.time)));
    elseif type_error == 2
        temp_2 = results1.x.moments(1,:)-((1-threshold)*ground_truth.x.moments(1,1:length(ground_truth.time)));
    end

    indexOfMinDiff = find(temp_2 < 0, 1);
    time_2_a = results1.time(indexOfMinDiff);

    % Check if there are disagreement times
    if (~isempty(time_1_a) && ~isempty(time_2_a)) 
        if time_1_a < time_2_a
            time1_mean = time_1_a;
        else
            time1_mean = time_2_a;
        end
    elseif (~isempty(time_1_a) && isempty(time_2_a))
        time1_mean = time_1_a;
    elseif (isempty(time_1_a) && ~isempty(time_2_a))
        time1_mean = time_2_a;
    else
        time1_mean = NaN;
    end

    %% Find the disagreement times for results1 for the variance
    if type_error == 1
        temp_1 = (ground_truth.x.moments(2,1:length(ground_truth.time))+threshold*ones(1,length(ground_truth.time)))-results1.x.moments(2,:);
    elseif type_error == 2
        temp_1 = (1+threshold)*ground_truth.x.moments(2,1:length(ground_truth.time))-results1.x.moments(2,:);
    end
    
    indexOfMinDiff = find(temp_1 < 0, 1);
    time_1_b = results1.time(indexOfMinDiff);

    if type_error == 1
        temp_2 = results1.x.moments(2,:)-(ground_truth.x.moments(2,1:length(ground_truth.time))-threshold*ones(1,length(ground_truth.time)));
    elseif type_error == 2
        temp_2 = results1.x.moments(2,:)-((1-threshold)*ground_truth.x.moments(2,1:length(ground_truth.time)));
    end

    indexOfMinDiff = find(temp_2 < 0, 1);
    time_2_b = results1.time(indexOfMinDiff);

    % Check if there are disagreement times
    if (~isempty(time_1_b) && ~isempty(time_2_b)) 
        if time_1_b < time_2_b
            time1_var = time_1_b;
        else
            time1_var = time_2_b;
        end
    elseif (~isempty(time_1_b) && isempty(time_2_b))
        time1_var = time_1_b;
    elseif (isempty(time_1_b) && ~isempty(time_2_b))
        time1_var = time_2_b;
    else
        time1_var = NaN;
    end

    %% Find the disagreement times for results2 for the mean
    if type_error == 1
        temp_1 = (ground_truth.x.moments(1,1:length(ground_truth.time))+threshold*ones(1,length(ground_truth.time)))-results2.x.moments(1,:);
    elseif type_error == 2
        temp_1 = (1+threshold)*ground_truth.x.moments(1,1:length(ground_truth.time))-results2.x.moments(1,:);
    end
    
    indexOfMinDiff = find(temp_1 < 0, 1);
    time_1_c = results2.time(indexOfMinDiff);

    if type_error == 1
        temp_2 = results2.x.moments(1,:)-(ground_truth.x.moments(1,1:length(ground_truth.time))-threshold*ones(1,length(ground_truth.time)));
    elseif type_error == 2
        temp_2 = results2.x.moments(1,:)-((1-threshold)*ground_truth.x.moments(1,1:length(ground_truth.time)));
    end

    indexOfMinDiff = find(temp_2 < 0, 1);
    time_2_c = results2.time(indexOfMinDiff);

    % Check if there are disagreement times
    if (~isempty(time_1_c) && ~isempty(time_2_c)) 
        if time_1_c < time_2_c
            time2_mean = time_1_c;
        else
            time2_mean = time_2_c;
        end
    elseif (~isempty(time_1_c) && isempty(time_2_c))
        time2_mean = time_1_c;
    elseif (isempty(time_1_c) && ~isempty(time_2_c))
        time2_mean = time_2_c;
    else
        time2_mean = NaN;
    end

    %% Find the disagreement times for results2 for the variance
    if type_error == 1
        temp_1 = (ground_truth.x.moments(2,1:length(ground_truth.time))+threshold*ones(1,length(ground_truth.time)))-results2.x.moments(2,:);
    elseif type_error == 2
        temp_1 = (1+threshold)*ground_truth.x.moments(2,1:length(ground_truth.time))-results2.x.moments(2,:);
    end
    
    indexOfMinDiff = find(temp_1 < 0, 1);
    time_1_d = results1.time(indexOfMinDiff);

    if type_error == 1
        temp_2 = results2.x.moments(2,:)-(ground_truth.x.moments(2,1:length(ground_truth.time))-threshold*ones(1,length(ground_truth.time)));
    elseif type_error == 2
        temp_2 = results2.x.moments(2,:)-((1-threshold)*ground_truth.x.moments(2,1:length(ground_truth.time)));
    end

    indexOfMinDiff = find(temp_2 < 0, 1);
    time_2_d = results1.time(indexOfMinDiff);

    % Check if there are disagreement times
    if (~isempty(time_1_d) && ~isempty(time_2_d)) 
        if time_1_d < time_2_d
            time2_var = time_1_d;
        else
            time2_var = time_2_d;
        end
    elseif (~isempty(time_1_d) && isempty(time_2_d))
        time2_var = time_1_d;
    elseif (isempty(time_1_d) && ~isempty(time_2_d))
        time2_var = time_2_d;
    else
        time2_var = NaN;
    end

end