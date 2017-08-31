%% Calculate slope and intersect for the initial region of each trace.
% Slope detection length is fixed to the value `init_slope_length' defined in the constants file.


function [intersct_array, coefs_array, slopes_array] = calculate_slope_and_intersect(ms2_combined, ms2_count, forced_start_nc_time)

%% Constants
constants;



% For each ms2 spots extract the data corresponding to this time, and fit it with a straight line.
intersct_array = zeros(1, ms2_count);
coefs_array = zeros(ms2_count, 2);
slopes_array = zeros(1, ms2_count);
% Parse collected traces
for i = 1: ms2_count
    cur_frames = ms2_combined(i).Frame;
	cur_fluo = ms2_combined(i).Fluo;
	frame_count = length(cur_frames);
    
    %% Detect the start of the nc
	% Cut fluorescence before the expected start of the nc. The traces should have been already manually adjusted to lie to the right of this value.
    cur_fluo (cur_frames < forced_start_nc_time/mins_per_frame) = 0;
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
        intersct_array(i) = intersct;
        slopes_array (i) = coefs(1)/fluo_per_polymerase;    % These slopes will be in pols/min
    end;
end;