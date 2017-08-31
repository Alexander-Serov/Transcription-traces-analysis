%% Calculate slope and intersect for the initial region of each trace.
% Slope detection length is fixed to the value `init_slope_length' defined in the constants file.
%
% The parameter `bl_detect_nc_start' defines whether the code will look for the beginning of the nuclear cycle. If it is 0, the slope will be measured in the
% first several minutes of the pre-defined nuclear cycle, as defined in the constants file.


function [intersct_array, coefs_array, slopes_array] = calculate_slope_and_intersect(ms2_data, ms2_count, forced_start_nc_time, mins_per_frame, bl_detect_nc_start)

%% Constants
constants;



% For each ms2 spots extract the data corresponding to this time, and fit it with a straight line.
intersct_array = zeros(1, ms2_count);
coefs_array = zeros(ms2_count, 2);
slopes_array = zeros(1, ms2_count);
% Parse collected traces
for i = 1: ms2_count
    cur_frames = ms2_data(i).Frame;
	cur_fluo = ms2_data(i).Fluo;
	frame_count = length(cur_frames);
% 	last_frame = max(cur_frames);
	
	% Cut fluorescence before the expected start of the nc. The traces should have been already manually adjusted to lie to the right of this value.
	cur_fluo (cur_frames < forced_start_nc_time/mins_per_frame) = 0;
    
	if bl_detect_nc_start
		%% Detect the start of the nc
		
		% Find the first non-zero frame
		slope_start_index = find(cur_fluo > 0, 1);
		% Skip trace if no fluorescence recorded in the nc
		if isempty(slope_start_index), continue, end;

		%% Detect the end of the slope
		% Assuming that ncs are correctly labeled.
		% Fixed slope length is used
		slope_end_index = slope_start_index + ceil(init_slope_length / mins_per_frame);
		% Choose the smallest of the pre-defined slope length or trace length
		slope_end_index = min (slope_end_index, frame_count);
	else
		cur_time = ms2_data(i).Time;
		% Find the frame corresponding to the pre-defined beginning of the nuclear cycle
		slope_start_index_1 = find(cur_time >= forced_start_nc_time, 1);
		
		% Make sure to start with non-zero data
		slope_start_index_2 = find(cur_fluo > 0, 1);
		% Skip trace if no fluorescence recorded in the nc or the beginning of the nuclear cycle is not included
		if isempty(slope_start_index_1) || isempty(slope_start_index_2)
			continue;
		else
			% Chose the later one to avoid empty data
			slope_start_index = max(slope_start_index_1, slope_start_index_2);
		end;
						
		% Set the end of the slope to a fixed interval after the start of the nuclear cycle
		slope_end_index = find(cur_time > forced_start_nc_time + init_slope_length, 1) - 1;	
		% Make sure that it's shorter than the trace length
		if isempty(slope_end_index)
			slope_end_index = frame_count;
		end;
		
		if slope_end_index == 0
			1;
		end;
		disp([slope_start_index, slope_end_index]);
	end;
	
	
	
	

    %% Fit the slope with a straight line.
    slope_frames_number = slope_end_index - slope_start_index + 1;
    
    if slope_frames_number > 1
		% Extract slope points.
		fluo_points_slope = cur_fluo(slope_start_index : slope_end_index);
		time_points_slope = cur_frames(slope_start_index : slope_end_index) * mins_per_frame;
	
        % Fit with a straight line.
		coefs = polyfit(time_points_slope, fluo_points_slope, 1);
		% Collect coefficients for all slopes, in fluo/min.
        coefs_array (i, :) = coefs;
        % Calculate intersect to the start of the nuclear cycle.
        intersct = - coefs(2) / coefs(1);
        intersct_array(i) = intersct;	% is probably in mins
        slopes_array (i) = coefs(1)/fluo_per_polymerase;    % These slopes will be in pols/min
    end;
end;



