function [neuronal_activity,figID] = HR_neuronal_activity(potential,time,x_R,figID)
% This function computes the spike count, duty cycle and ISI of the time
% series potential having resting potential x_R (not updated, transition
% time still considered)

    % Temporary variables for successive resizing along iterations, the
    % part of the signal already analysed is discarded
    if size(potential,2) > size(potential,1)
        potential_temp = potential';
    else
        potential_temp = potential;
    end

    % potential_temp = smooth(potential_temp);   % default span is 5

    if size(time,2) > size(time,1)
        time_temp = time';
    else
        time_temp = time;
    end

    % Initialisation of iteration variables
    burst_count = 0;    % number of full bursts in the given signal
    %stop_flag = 0;      % indicates when to stop the cycle
    first_time = 1;
    spiking = 0;        % indicates when there is just spiking
    time_indx = 0;      % tracks the progression along time

    glob_min = min(potential_temp);   % for discarding parts of the signal
    
    % Initialisation of performance indices
    spike_count = [];
    ISI = [];
    max_amplitude = 0;
    t_hyperpol = [];
    duty_cycle = [];
    t_period = [];

    % Arbitrary settings
    prominence = 0.07;
    rescaling = 0.25;   % originally 0.25

    while 1
        % Check first if the signal is below the resting potential, if it is then skip the first part of the signal
        if potential_temp(1) < x_R && first_time
            first_time = 0;
            begin_next_burst_indx = find(potential_temp >= x_R, 1); 
            potential_temp = potential_temp(begin_next_burst_indx:end);
            time_temp = time_temp(begin_next_burst_indx:end);
            time_indx = time_indx + begin_next_burst_indx;
            continue;
        end

        % Computation of the first index at which the potential hyperpolarizes after the burst offsert
        end_burst_indx = find(potential_temp <= x_R, 1);
        
        %--------------------------------------------------------------------------------------------------------------------------------------------------
        %--------------------------------------------------------------------------------------------------------------------------------------------------
        if isempty(end_burst_indx) % either the time sequence is finished or there is no bursting after some time

            if burst_count ~= 0 && time_indx > length(time)/2
                if spiking
                    disp(['The time sequence contains ',num2str(burst_count),' full spiking cycles.']);
                    break;
                else
                    disp(['The time sequence contains ',num2str(burst_count),' full bursting cycles.']);
                    break;
                end
            else
                t_hyperpol = 0;    % there is no hyperpolarization in the regime
                spiking = 1;
                % rescale the resting potential since it is never crossed and restrt the iteration
                x_R = glob_min + rescaling*abs(glob_min);    % this value will be crossed for sure, so the stop_flag will be reached eventually
                disp('There is no bursting activity, the potential is always above the resting level.');
                continue
            end
            
            % After this last iteration it is possible to leave the cycle
            %stop_flag = 1;
        end

        % % Exit the cycle only when the signal analysis is completed
        % if stop_flag
        %     break;
        % end
        
        %-----------------------------------------------------------------------------------------------------------------------------------------------------
        % -------------------------- Split the signal in bursting periods and modify it to count the peaks ---------------------------------------------------
        %-----------------------------------------------------------------------------------------------------------------------------------------------------
        

        %--------------------------------------------------------------------------------------------------------------------------------------------------------
        % Check if additional iterations are necessary
        % if isempty(end_burst_indx)
        %     burst = potential_temp;                            % burst is all that is left of the signal
        % else
        %---------------------------------------------------------------------------------------------------------------------------------------------------------

        burst = potential_temp(1:end_burst_indx);          % burst is a certain portion of the signal
        time_indx = time_indx + end_burst_indx;            % time index is updated
        

        %----------------------------------------------------------------------------------------------------------------------------------------------------------
        % end
        %----------------------------------------------------------------------------------------------------------------------------------------------------------


        % Remove the curve from the resting potential to the first peak in order to avoid having local maxima in the signal

        % Find first the index of the first peak
        peaks_vec_original = islocalmax(burst,'MinProminence',prominence);          % positions of local maxima in the original signal piece (default prominence is zero and default minimum separation is 0)
        max_peaks_indx = find(peaks_vec_original);

        if isempty(max_peaks_indx)
            potential_temp = potential_temp(end_burst_indx:end);
            time_temp = time_temp(end_burst_indx:end);
            time_indx = time_indx + end_burst_indx;
            begin_next_burst_indx = find(potential_temp >= x_R, 1); 
            potential_temp = potential_temp(begin_next_burst_indx:end);
            time_temp = time_temp(begin_next_burst_indx:end);
            time_indx = time_indx + begin_next_burst_indx;
            continue;
        end

        max_peaks = maxk(burst(peaks_vec_original),sum(peaks_vec_original));   % vector containing the values of the local maxima
        for i=1:length(max_peaks_indx)                                         % first peak is within the 4 largest peaks
            if burst(max_peaks_indx(i)) == max_peaks(1) || ...
               burst(max_peaks_indx(i)) == max_peaks(2) || ...
               burst(max_peaks_indx(i)) == max_peaks(3) || ...
               burst(max_peaks_indx(i)) == max_peaks(4)
                first_peak = max_peaks_indx(i);
                break;
            end
        end
        burst(1:first_peak-1) = glob_min;

        max_amplitude = max([max_amplitude,max(max_peaks)]);   % compute max amplitude of all the spikes

        % Smooth the curve at the burst offset until the signal reaches the
        % resting potential
        n_valleys = islocalmin(burst,'MinProminence',prominence);
        n_valleys(1:first_peak) = 0;                           % do not considered the part of the signal before the first spike
        n_valleys(end) = 0;                                    % always discard the very last value at the resting potential
        valleys_indx = find(n_valleys);
        min_valley = min(burst(valleys_indx));                 % minimum value of the potential during bursting
        
        if ~isempty(min_valley)     % typically it is empty for tonic spiking behaviour
            burst(burst<min_valley) = glob_min;
        end

        % Compute the local maxima locations in the smoothed vector and count the number of spikes
        peaks_vec = islocalmax(burst,'MinProminence',prominence);

        % if isempty(find(burst==glob_min,1))                    % discard the maximum before the discontinuity in burst (if any)
        %     peaks_vec(find(burst==glob_min,1)-1) = 0;  
        % end

        spike_count = [spike_count,sum(peaks_vec)];
        peaks_indx = find(peaks_vec);

        if isempty(peaks_indx)                                 % we are at the end of the signal and in that portion there are no local maxima
            break;
        end
        
        last_peak = peaks_indx(end);                           % index of last spike for computation of the duty cycle
        
        %---------------------------------------------------------------------------------------------------------------------------------------------------
        % ------------------------------------------ Compute the times of the signal portion ---------------------------------------------------------------
        %---------------------------------------------------------------------------------------------------------------------------------------------------
        % Compute the time of bursting for the duty cycle
        t_bursting = time_temp(last_peak)-time_temp(first_peak);    

        % Compute the ISI
        ISI_temp = diff(time_temp(peaks_indx));


        %-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % if burst_count == 0
        %     ISI = [ISI_temp(1:end)];
        %     % Add the time period from the last spike until the potential
        %     % reaches the resting potential
        %     ISI = [ISI;time_temp(end_burst_indx)-time_temp(last_peak)];
        %     if isempty(end_burst_indx)
        %         t_period_i = sum(ISI_temp)+time_temp(end)-time_temp(peaks_indx(end));   % resting potential is not reached, we add the residual time of the signal
        %     else
        %         t_period_i = sum(ISI_temp)+time_temp(end_burst_indx)-time_temp(last_peak);
            % end
        % else
            % Add to the last ISI the bit from the potential reaching x_R
            % again after hyperpolarization to the first spike of the new
            % bursting sequence
            %ISI(end) = ISI(end)+(time_temp(peaks_indx(1))-time_temp(1));
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


        % Add the time period from the last spike until the potential
        % reaches the resting potential
        ISI = [ISI;ISI_temp;time_temp(end_burst_indx)-time_temp(peaks_indx(end))];


        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        % if isempty(end_burst_indx)
        %     t_period_i = (time_temp(peaks_indx(1))-time_temp(1))+sum(ISI_temp)+(time_temp(end)-time_temp(last_peak));  % resting potential is not reached, we add the residual time of the signal
        % else
        %-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        t_period_i = (time_temp(peaks_indx(1))-time_temp(1))+sum(ISI_temp)+(time_temp(end_burst_indx)-time_temp(last_peak));
        k_temp = time_temp(peaks_indx(1))-time_temp(1);

        %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
        %end
        % end
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        %------------------------------------------------------------------------------------------------------------------------------------------------------
        %-------------------------------------------- Resize the potential time series for the next iteration -------------------------------------------------
        %------------------------------------------------------------------------------------------------------------------------------------------------------
        potential_temp = potential_temp(end_burst_indx:end);
        time_temp = time_temp(end_burst_indx:end);

        begin_next_burst_indx = find(potential_temp >= x_R, 1);  % start next iteration when the potential rises from the resting potential
        
        % The system converges to a stable node, typical for low values of b
        if isempty(begin_next_burst_indx)
            if spiking
                disp(['The time sequence contains ',num2str(burst_count+1),' full spiking cycles.']);
            else
                t_hyperpol = [t_hyperpol,Inf];  % if there is spiking the hyperpolarization time is already set to zero
                disp(['The time sequence contains ',num2str(burst_count+1),' full bursting cycles.']);
            end
            disp('After the last burst the system remains below the resting level for the rest of the time series.');
            ISI(end) = ISI(end) + (time(end)-time_temp(1));
            break;
        end

        time_indx = time_indx + begin_next_burst_indx;
        
        % Add to the last ISI the hyperpolarization period 
        t_hyperpolarized = time_temp(begin_next_burst_indx)-time_temp(1);
        t_period_i = t_period_i + t_hyperpolarized;
        if ~spiking
            t_hyperpol = [t_hyperpol,t_hyperpolarized];             
        end
        ISI(end) = ISI(end) + t_hyperpolarized;

        % Compute the duty cycle before resizing time and potential
        duty_cycle = [duty_cycle,t_bursting/t_period_i];

        % Remove the hyperpolarization period for the next iteration
        potential_temp = potential_temp(begin_next_burst_indx:end);
        time_temp = time_temp(begin_next_burst_indx:end);
        
        % A full burst cycle is completed
        burst_count = burst_count+1;
        t_period = [t_period,t_period_i];

        % Do not count the hyperpolarization period
        if t_bursting == 0
            ISI(end) = ISI(end)+k_temp;
        else
            ISI = ISI(1:end-1);
        end

        if first_time   % discard the very first part of the signal
            first_time = 0;
            spike_count = [];
            ISI = [];
            max_amplitude = 0;
            t_hyperpol = [];
            duty_cycle = [];
            t_period = [];
            burst_count = 0;
        end

    end

    % Plot of neuronal activity
    figID = figID+1;
    figure(figID);
    histogram(ISI,'BinMethod','sqrt');
    title('ISI histogram');

    figID = figID+1;
    figure(figID);
    histogram(spike_count,'BinMethod','sqrt');
    title('Spike count histogram');

    figID = figID+1;
    figure(figID);
    histogram(duty_cycle,'BinMethod','sqrt');
    title('Duty cycle histogram');
    
    % Output struct
    neuronal_activity.spike_count = spike_count;
    neuronal_activity.ISI = ISI;
    neuronal_activity.burst_count = burst_count;
    neuronal_activity.duty_cycle = duty_cycle;
    neuronal_activity.max_amplitude = max_amplitude;
    neuronal_activity.t_period = t_period;
    neuronal_activity.t_hyperpol = t_hyperpol;


end

