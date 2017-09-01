

function analysis_7_new_data(gene_ind, dataset_ind, nuc_cyc)


%% Constants
constants;



%% Initialize
% Calculate theoretical slope
theoretical_slope = k/(1+sqrt(l))^2;	% k in bp/min



% % %% Cycling through all the data
% % % I am not sure why I reduced the cycle

% gene_ind = 1;
%% Load gene and construct name
gene_name = gene_names_array{gene_ind};
% % dataset_ind = 1;
dataset_name = dataset_names_array{dataset_ind};
fprintf('Analysis in progress. Gene: %i/%i. Construct: %i/%i. Nuclear cycle: %i.\n', gene_ind, 3, dataset_ind, 3, nuc_cyc);
        
% % Skip Knirps, no shadow in nc13 because no data
% if gene_ind == 2 && dataset_ind == 3 && nuc_cyc == 13, continue, end;



%% Defining different bin widths for different genes
if strcmp(gene_name, 'HunchBack') || strcmp(gene_name, 'SNAIL')
    bin_width = 0.05;
elseif strcmp(gene_name, 'Knirps')
    bin_width = 0.05;
end;



%% Defining nuc. cyc.-specific parameters
if nuc_cyc == 13
    forced_start_nc_time = forced_start_nc_time_13;  % in mins
    xlim_vector = [32, 55];
elseif nuc_cyc == 14
    forced_start_nc_time = forced_start_nc_time_14;  % in mins
    xlim_vector = [47, 60];
end;



%% Load and align data
[Data, embryos, dataset_count, time_alignment_mins] = choose_datasets_and_time_alignment(gene_name, dataset_name, nuc_cyc);



%% Combine, pre-align and filter data from multiple embryos for the given gene, construct and nuclear cycle

% Initialize
ms2_count = 0;
ms2_input_count = 0;
mins_per_frame = orig_mins_per_frame;

% Parse all available embryos (data sets)
for embryo_number = embryos
    % Parse through all traces recorded in a given embryo
    for i = 1: length(Data(embryo_number).ms2)
		ms2_input_count = ms2_input_count + 1;
        % Get current trace
        cur_trace = Data(embryo_number).ms2(i);
        
        % Check that the nuclear cycle is correctly recorded
        if cur_trace.nuc_cyc ~= nuc_cyc, continue, end;
           
        % Skip if the trace has only one data point
        if length(cur_trace.Frame) <= 1, continue, end;
		
		% Skip trace, if some of the AP positions are negative! (Who prepared these data??)
		if min(cur_trace.APpos) < 0, continue, end;
        
        % Calculate trace integral
        observation_integral = trapz(cur_trace.Frame(:), cur_trace.Fluo(:));
		% Convert to mins and polymerase number
		observation_integral = observation_integral * mins_per_frame / fluo_per_polymerase;
        % Threshold based on the integral fluorescence value or if we see negative fluorescence anywhere in the trace
        if observation_integral <= integral_threshold_value || min(cur_trace.Fluo(:)) < 0, continue, end;
        
        % Make a rough (manual) adjustment of the time
        cur_trace.('Frame') = cur_trace.('Frame') + round(time_alignment_mins(embryo_number) / mins_per_frame);
        
        % Increase the counter and save trace to the common array
        ms2_count = ms2_count + 1;
		ms2_combined(ms2_count) = cur_trace;    
    end;
end
fprintf('Kept %i filtered traces out of %i for %s %s in nc %i.\n', ms2_count, ms2_input_count, gene_name, dataset_name, nuc_cyc);


% If no traces found, skip case
if ms2_count == 0
	fprintf('No traces available for %s %s in nc %i. Skipping case.\n', gene_name, dataset_name, nuc_cyc);
	return;
end;

    

%% Plot roughly synced data for the current gene, construct and nuclear cycle

figure(1);
clf;
hold on;
for i = 1:ms2_count
    plot (ms2_combined(i).Frame * mins_per_frame, ms2_combined(i).Fluo/fluo_per_polymerase);
end;

% Label plot
xlabel('Time, min');
ylabel('Number of active polymerases');
title('Roughly adjusted fluorescence data');
hold off;
pause(0.1);
    


%% Calculate slope and intersect for the initial region of each trace.
% The slope is calculated on the first several frames starting from the first non-zero point.
% The number of frames to consider is pre-defined in the constants file, and is the same for all traces

[intersct_array, slopes_array, slope_start_frame_array] = calculate_slope_and_intersect_smart(ms2_combined, ms2_count, mins_per_frame);



%% Filter data based on slope and intersect, then create a refined mesh and apply time shifts when projecting data.
% The unifrom mesh is needed for the calculation of averages.
% There will be nothing before the start of the cycle, because we synchronize everything by the first few frames.
% This may creast problem if data contain noise before nc. We try to filter that by ignoring negative slopes. Very low slopes should also be visible later.

% Calculate max. frame number (maximum time) for the collected traces
max_frame_number = max(ms2_combined(1).Frame);
for i = 1:ms2_count
    max_frame_number = max(max_frame_number, max(ms2_combined(i).Frame));
end;
min_frame_number = 1;

% Create refined mesh.
original_time_mesh = (1:max_frame_number) * mins_per_frame;
time_step = desired_time_step;    % in mins
% Extend the time frame by the allowed shift
time_mesh = new_time_start: time_step : new_time_end;
time_mesh_length = length(time_mesh);

% Initialize new data structure to store filtered traces.
new_data = struct();
new_data_count = 0;

% Apply the calculated shifts, while keeping only traces with positive slopes.
for i=1:ms2_count
    current_shift = forced_start_nc_time - intersct_array(i);
    if slopes_array(i) <= 0
        slopes_array(i) = NaN;
        continue;
    end;
    
    % Extract data to convert.
	old_time = ms2_combined(i).Frame * mins_per_frame;  % in mins
    old_fluo = ms2_combined(i).Fluo;	% in fluo
    
	% Apply a shift
    old_time = old_time + current_shift;
   
    % Determine new frame numbers on which to project the data points
	% Find the new starting frame
	% Only use interpolation
    new_start_frame = find(time_mesh >= old_time(1), 1);
	new_start_time = time_mesh(new_start_frame);
    
	% Find the new end frame
    new_end_frame = find(time_mesh > old_time(end), 1) - 1;
	new_end_time = time_mesh(new_end_frame);
    
    new_trace_time_points = new_start_time : time_step : new_end_time;
	
    % Get linearly interpolated fluorescence data on the new time points
    fluo_at_new_trace_time_points = interp1(old_time, old_fluo, new_trace_time_points, 'linear');
    
    % Save results to new data structure.
	new_data_count = new_data_count + 1;
    new_data(new_data_count).Frame = new_start_frame:new_end_frame;
    new_data(new_data_count).Time = new_trace_time_points;
    new_data(new_data_count).Fluo = fluo_at_new_trace_time_points;
	new_data(new_data_count).AP = ms2_combined(i).APpos;
	new_data(new_data_count).mean_AP = mean(ms2_combined(i).APpos);
	new_data(new_data_count).y_pos = ms2_combined(i).ypos;
	new_data(new_data_count).mean_y_pos = mean(ms2_combined(i).ypos);
	new_data(new_data_count).Slope = slopes_array(i);
	
	% Save slope normalized to sMax
    new_data(new_data_count).norm_slope = slopes_array(i) / theoretical_slope;
end;

% Update the frame rate.
mins_per_frame = time_step;



%% Plot traces accurately shifted.
% The beginning of the ncs should have been set already to exactly the value defined in the constants file.

figure(2);
clf;
hold on;
for i = 1:new_data_count
	% Get data.
	cur_time = new_data(i).Time;
	cur_fluo = new_data(i).Fluo / fluo_per_polymerase;
	% Plot.
	plot (cur_time, cur_fluo);
end;

% Adjust limits.
xlim(xlim_vector);
ylim([0,110]);
% Label.
xlabel('Time, min');
ylabel('Number of active polymerases');
title('Synced fluorescence data');

% Print statistics
fprintf('Slope identified in %i traces out of %i original traces for %s %s in nc %i.\n', new_data_count, ms2_input_count, gene_name, dataset_name, nuc_cyc);

% Adding the theoretical prediction
% hold on;
theor_time_mesh = forced_start_nc_time:0.05:(forced_start_nc_time + 4);
theor_slope_func = @(t) k/(1+sqrt(l))^2 * t;		% k in min^{-1}, t in min
theor_slope_data = theor_slope_func(theor_time_mesh - forced_start_nc_time);
plot(theor_time_mesh, theor_slope_data, '--', 'LineWidth', 2, 'color', 'black');

% Adding the maximal particle number
temp_time_mesh = xlim_vector(1):0.05:xlim_vector(2);
plot(temp_time_mesh, ones(1,length(temp_time_mesh)) * L/(l+sqrt(l)), ':', 'color', 'black',...
    'LineWidth', 2);
hold off;
pause(0.1);



%% Group remaining traces into AP bins

% Initialize
ms2_mean_AP = zeros(1, new_data_count);
new_data_to_bin_table = zeros(1, new_data_count);

% Get the average AP position of each trace
for i=1:new_data_count
    ms2_mean_AP(i) = new_data(i).mean_AP;
end;

% If the AP limits are not manually defined, use the min and max of the observed values
if ~manual_bins
    min_ms2_AP = min(ms2_mean_AP);
    max_ms2_AP = max(ms2_mean_AP);
end;

% Define bin locations
bins_borders = (min_ms2_AP - bin_width / 2) : bin_width : (max_ms2_AP + bin_width * 3 / 2);
bins_count = length(bins_borders) - 1;
bins_centers = (bins_borders(1 : end - 1) + bins_borders(2 : end)) / 2;

% Assign each trace to a bin
for  i = 1:new_data_count
    current_bin = ceil((new_data(i).mean_AP - bins_borders(1)) / bin_width);
    % Save bin to table
	new_data_to_bin_table (i) = current_bin;
	% Save bin to trace
	new_data(i).bin = current_bin;
end;



%% Collect traces of each bin
% This 3D array is not very economic, but I don't have memory concerns for the moment.

bin_trace_count = zeros(bins_count, 1);
bin_fluo_data_new = ones(bins_count, time_mesh_length, new_data_count) * NaN;
bin_norm_slope_data_new = ones(bins_count, new_data_count) * NaN;

% Parse traces
for i = 1 : new_data_count
    bin = new_data(i).bin;
	% Get frames corresponding to the trace.
	frames = new_data(i).Frame;
	
	% Save fluo data to the new structure
	bin_trace_count(bin) = bin_trace_count(bin) + 1;
	bin_fluo_data_new(bin, frames, bin_trace_count(bin)) = new_data(i).Fluo(:);
	
	% Save normalized slope
	bin_norm_slope_data_new(bin, bin_trace_count(bin)) = new_data(i).norm_slope;
end;



%% Calculating mean and STD across traces for each frame and bin

bin_fluo_data_mean = mean(bin_fluo_data_new, 3, 'omitnan');
bin_fluo_data_STD = std(bin_fluo_data_new, 0, 3, 'omitnan');
% Count the number of traces per bin and frame
bin_frame_trace_count = sum(~isnan(bin_fluo_data_new), 3);

% Calculate mean normalized slope
bin_norm_slope_mean = mean(bin_norm_slope_data_new, 2, 'omitnan');



%% Calculate bootstrap confidence intervals for fluo value and normalized slope.
% For each bin and frame bootstrap to find confidence intervals.
% Some frames will be skipped to accelerate calculations.

disp('Calculating bootstrap confidence intervals');

% Initialize
bin_fluo_data_shifted_confidence_intervals = ones(bins_count, time_mesh_length, 2) * NaN;
bin_norm_slope_confidence_intervals = ones(bins_count, 2) * NaN;
tic;
% Define the bootstrapping function
bootfunc = @(x) mean(x);
bs_samples_count = bootstrap_samples_count;

for bin = 1 : bins_count
    %% Bootstrap fluorescence in each frame
	parfor frame = 1:time_mesh_length
		
		% Extract data
		cur_fluo = bin_fluo_data_new(bin, frame, :);
		% Remove NaN values
		cur_fluo = cur_fluo(~isnan(cur_fluo));
		
		% Skip if has less than 2 elements
		if length(cur_fluo) < 2, continue, end;
   
		% Calculate.
		cur_conf_int = bootci(bs_samples_count, bootfunc, cur_fluo);

		%% Convert absolute confidence intervals to relative intervals
% 		% Keep the left boundary non-negative
% 		if bin_fluo_data_shifted_confidence_intervals(bin, frame, 1) > 0
		cur_conf_int(1) = bin_fluo_data_mean(bin, frame) - cur_conf_int(1);
% 		else
% 			bin_fluo_data_shifted_confidence_intervals(bin, frame, 1) = 0;
% 		end;
		
% 		if bin_fluo_data_shifted_confidence_intervals(bin, frame, 2)>0
		cur_conf_int(2) = cur_conf_int(2) - bin_fluo_data_mean(bin, frame);
% 		else
% 			bin_fluo_data_shifted_confidence_intervals(bin, frame, 2) = 0;
% 		end;
		bin_fluo_data_shifted_confidence_intervals(bin, frame, :) = cur_conf_int;
	end;
    
    %% Bootstrapping the normalized slopes in each bin
	% Extract data
	cur_slopes = bin_norm_slope_data_new(bin, :);
	% Remove NaN values
	cur_slopes = cur_slopes(~isnan(cur_slopes));
	
	% Skip if has less than 2 elements
	if length(cur_slopes) < 2, continue, end;
	
	% Calculate.
    bin_norm_slope_confidence_intervals(bin,:) = bootci(bs_samples_count, bootfunc, cur_slopes);
    
    fprintf('Bootstrapping progress: %.1f%%\n', bin / bins_count * 100);
end;

fprintf('Bootstrap confidence intervals calculation completed in %.2f s\n', toc);



%% Old plot functions currently unused.
% % % plot_old;


%% Output the considered AP interval
fprintf('AP interval for %s %s: [%.2f, %.2f]\n', gene_name, dataset_name, min_ms2_AP, max_ms2_AP);


% % % %% Let's find the maximal slope.
% % % % For this let's first find all slopes. 
% % % 
% % % % plot(new_time_mesh,...
% % % %         bin_fluo_data_shifted_mean(:,bin)/fluo_per_polymerase, 'LineWidth', 2);
% % % 
% % % bin_mean_slopes = zeros(1, bins_count);
% % % 
% % % % Selecting the right time interval
% % % time_interval_slope = 2;    % where to look for the slope, in mins
% % % start_time_mins = forced_start_nc_time;
% % % t_ind_slope_start = find(time_mesh > start_time_mins, 1, 'first');
% % % t_ind_slope_end = find(time_mesh > start_time_mins+time_interval_slope, 1, 'first');
% % % 
% % % indices_to_consider = t_ind_slope_start:t_ind_slope_end;
% % % 
% % % for bin = 1:bins_count
% % %     
% % %     % From the indices to consider I have to choose only those that are not
% % %     % NaN
% % %     indices_not_nan = find(~isnan(bin_fluo_data_shifted_mean(:,bin)));
% % %         
% % %     indices_to_consider_not_nan = intersect(indices_to_consider, indices_not_nan);
% % %     
% % %     
% % %     fit_result= ((time_mesh(indices_to_consider_not_nan)' - start_time_mins)\...
% % %         [ones(length(indices_to_consider_not_nan),1), bin_fluo_data_shifted_mean(indices_to_consider_not_nan,bin)]);
% % %     bin_mean_slopes(bin) = fit_result(2);
% % % end;
% % % 
% % % bin_mean_slopes = bin_mean_slopes/fluo_per_polymerase;



%% Print out the maximal slope over all bins
max_norm_experimental_slope = max(max(bin_norm_slope_data_new));
% theoretical_slope = k/(1+sqrt(l))^2;
fprintf('\nMaximal experimental slope: %.2f pol/min\n', max_norm_experimental_slope * theoretical_slope);
fprintf('Theoretical slope: %.2f pol/min\n', theoretical_slope);
fprintf('Percent of the theoretical slope: <= %.0f%%\n', max_norm_experimental_slope * 100);


%% Prepare data for output
slopes_array = ones(new_data_count, 1) * NaN;
norm_slopes_array = ones(new_data_count, 1) * NaN;
for trace = 1 : new_data_count
	slopes_array(trace) = new_data(trace).Slope;
	norm_slopes_array(trace) = new_data(trace).norm_slope;
end;


%% Output the processed data
output_folder = './processed_data/';
output_filename = sprintf('%s_%s_nc_%i.mat', gene_name, dataset_name, nuc_cyc);
output_full_path = strcat(output_folder, output_filename);

save(output_full_path, 'bin_fluo_data_mean', 'time_mesh', 'time_step', 'time_mesh_length', 'bins_borders', 'bins_count', 'bins_centers', 'bin_width', ...
        'bin_fluo_data_shifted_confidence_intervals', 'bin_fluo_data_STD', 'bin_fluo_data_new', 'bin_frame_trace_count', 'bin_trace_count',...
        'slopes_array', 'norm_slopes_array', 'bin_norm_slope_data_new', 'bin_norm_slope_mean', 'bin_norm_slope_confidence_intervals', 'new_data',...
        'new_data_count', 'forced_start_nc_time', 'intersct_array', 'ms2_input_count');










