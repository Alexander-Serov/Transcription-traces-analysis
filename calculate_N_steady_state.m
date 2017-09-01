


function steady_state_number = calculate_N_steady_state(ms2_combined, slopes_array, mins_per_frame, forced_start_nc_time,...
    intersct_array, init_slope_length, half_width_max_rgn_ind, fluo_per_polymerase)


% This code is modified to give a better estimate of the plateau value.
% Now we will take not just the peak value, but also 2 points on each
% side of it and calculate the average
ms2_count = size(ms2_combined, 2);
steady_state_number = zeros(ms2_count, 1);
for ms2_ind=1:ms2_count
    % Taking only traces for which we have a well defined slope
    if ~isnan(slopes_array(ms2_ind))
        % Finding index max. value
        cur_trace_length = length(ms2_combined(ms2_ind).Fluo);
        max_ind = find(ms2_combined(ms2_ind).Fluo == max(ms2_combined(ms2_ind).Fluo));
        % Determining plateau start index
        plateau_rgn_frame_indices = ((ms2_combined(ms2_ind).Frame * mins_per_frame) + (forced_start_nc_time - intersct_array(ms2_ind))) - ...
            (forced_start_nc_time + init_slope_length) > 0;
        plateau_start_index = find(plateau_rgn_frame_indices, 1);
        % Chossing a maximum of 5 indices around the max value to include into the mean plateau
        % calculation
        if plateau_start_index > 0
            % Using the plateau value only if the plateau region was detected
            rgn_start_index = max([max_ind - half_width_max_rgn_ind, plateau_start_index]);
            rgn_end_index = min([max_ind + half_width_max_rgn_ind, cur_trace_length]);
            % Calculating the plateau estimate
            cur_data = ms2_combined(ms2_ind).Fluo;
            cur_data = cur_data(rgn_start_index:rgn_end_index)/fluo_per_polymerase;
            steady_state_number(ms2_ind) = mean(cur_data);
        else
            % Ignoring data if the plateau region wasn't detected
            steady_state_number(ms2_ind) = NaN;
        end;
    else
        % If the slope is not defined, saving as NaN
        steady_state_number(ms2_ind) = NaN;
    end;
end;


