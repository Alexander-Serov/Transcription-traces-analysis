%% Calculate slope and intersect for the initial region of each trace.
% Slope detection length is fixed to the value `init_slope_length_frames' defined in the constants file.



function [intersct_array, slopes_array, slope_start_frame_array] = calculate_slope_and_intersect_smart(ms2_data, ms2_count, mins_per_frame)

%% Constants
constants;



% For each ms2 spots extract the data corresponding to this time
intersct_array = -1 * ones(1, ms2_count);
slopes_array = -1 * ones(1, ms2_count);
slope_start_frame_array = -1 * ones(1, ms2_count);

% Parse provided traces. The operation is fast, so no need to parallelize
for trace = 1: ms2_count
	% Extract current trace
	cur_frames = ms2_data(trace).Frame;
	cur_fluo = ms2_data(trace).Fluo;
	frame_count = length(cur_frames);
	
	% The initial slope length is fixed and pre-defined
	slope_length_frames = init_slope_length_frames;
	
	% Find the first data-containing frame
	start_frame = find(cur_fluo > 0, 1);
	% If no fluorescence found, skip trace
	if isempty(start_frame), continue, end;
		
	% If the expected slope goes outside recorded points, skip trace
	if start_frame + slope_length_frames > frame_count, continue, end;
	
	% Calculate a linear fit
	f = cur_fluo(start_frame : start_frame + slope_length_frames) / fluo_per_polymerase;
	t = cur_frames(start_frame : start_frame + slope_length_frames) * mins_per_frame;
	slope = (mean(f .* t) - mean(f) * mean(t)) / (mean(t.^2) - mean(t)^2);	% in pol/min
	b = mean(f) - slope * mean(t);	
	intersct = - b/slope;	% in mins
	
	% Save slope and intersect
	slopes_array(trace) = slope;
	intersct_array(trace) = intersct;
	slope_start_frame_array(trace) = start_frame;
end;



