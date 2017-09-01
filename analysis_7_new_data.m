
set (0, 'DefaultAxesFontSize', 20);
clear;



%% Constants
constants;



%% Initialize
% Calculate theoretical slope
theoretical_slope = k/(1+sqrt(l))^2;	% k in bp/min



%% Cycling through all the data
% I am not sure why I reduced the cycle
for gene_ind = 1:1 %3
    gene_name = gene_names_array{gene_ind};
    for dataset_ind = 1:1 %3
        dataset_name = dataset_name_array{dataset_ind};
        fprintf('Analysis in progress. Gene: %i/%i. Construct: %i/%i\n', gene_ind, 3, dataset_ind, 3);
        
% Skip Knirps, no shadow in nc13 because no data
if gene_ind == 2 && dataset_ind == 3 && nuc_cyc == 13, continue, end;



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
[Data, datasets, dataset_count, time_alignment_mins] = choose_datasets_and_time_alignment(gene_name, dataset_name, nuc_cyc);



%% Combine, pre-align and filter data from multiple embryos for the given gene, construct and nuclear cycle

% Initialize
ms2_count = 0;
ms2_input_count = 0;
mins_per_frame = orig_mins_per_frame;

% Parse all available embryos (data sets)
for dataset_number = datasets
    % Parse through all traces recorded in a given embryo
    for i = 1: length(Data(dataset_number).ms2)
		ms2_input_count = ms2_input_count + 1;
        % Get current trace
        cur_trace = Data(dataset_number).ms2(i);
        
        % Check that the nuclear cycle is correctly recorded
        if cur_trace.nuc_cyc ~= nuc_cyc, continue, end;
           
        % Skip if the trace has only one data point
        if length(cur_trace.Frame) <= 1, continue, end;
        
        % Calculate trace integral
        observation_integral = trapz(cur_trace.Frame(:), cur_trace.Fluo(:));
		% Convert to mins and polymerase number
		observation_integral = observation_integral * mins_per_frame / fluo_per_polymerase;
        % Threshold based on the integral fluorescence value or if we see negative fluorescence anywhere in the trace
        if observation_integral <= integral_threshold_value || min(cur_trace.Fluo(:)) < 0, continue, end;
        
        % Make a rough (manual) adjustment of the time
        cur_trace.('Frame') = cur_trace.('Frame') + round(time_alignment_mins(dataset_number) / mins_per_frame);
        
        % Increase the counter and save trace to the common array
        ms2_count = ms2_count + 1;
		ms2_combined(ms2_count) = cur_trace;    
    end;
end
fprintf('Kept %i filtered traces out of %i for %s %s in nc %i.\n', ms2_count, ms2_input_count, gene_name, dataset_name, nuc_cyc);

    

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
time_mesh = (original_time_mesh(1) - time_shift_limit): time_step : (original_time_mesh(end) + time_shift_limit);
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
	new_data(new_data_count).Slope = slopes_array(i);
	
	% Save slope normalized to sMax
    new_data(new_data_count).Norm_Slope = slopes_array(i) / theoretical_slope;
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



%% Group remaining traces into AP bins

% Initialize
ms2_mean_AP = zeros(1, new_data_count);
new_data_to_bin_table = zeros(1, new_data_count);

% Get the average AP position of each trace
for i=1:new_data_count
    ms2_mean_AP(i) = mean(ms2_combined(i).APpos);
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
for  i =1:new_data_count
    current_bin = ceil((ms2_mean_AP(i) - bins_borders(1)) / bin_width);
    new_data_to_bin_table (i) = current_bin;
end;



%% Now plotting non-shifted and shifted data grouped in bins

% Initializing
% Non-shifted data
% frame_bins_count = ceil(max_frame_number-min_frame_number+1);

bin_fluo_data_non_shifted = cell(max_frame_number, bins_count);
% bin_fluo_data_count_non_shifted = cell(max_frame_number, bins_count);
% Shifted data
bin_fluo_data_shifted = cell(time_mesh_length, bins_count);
bin_normalized_slopes = cell(1, bins_count);
% bin_fluo_data_count_shifted = cell(new_time_mesh_length, bins_count);


% Collecting data in bins
for i = 1: ms2_count
    current_bin = new_data_to_bin_table(i);
    if current_bin>bins_count || current_bin<=0
        continue;
    end;
    
    % Non-shifted data
    for j=1:length(ms2_combined(i).Frame)
        current_frame = ms2_combined(i).Frame(j);
        if current_frame <=0
            continue;
        end;

        
        bin_fluo_data_non_shifted{current_frame, current_bin} = ...
            horzcat(bin_fluo_data_non_shifted{current_frame, current_bin}, ms2_combined(i).Fluo(j));
        
    end;
    
    % Shifted data
    % Checking if we have shifted data for this ms2
    if ~isempty(new_data.ms2{i})
        for j=1:length(new_data.ms2{i}.newFrame)
            current_frame = new_data.ms2{i}.newFrame(j);

            bin_fluo_data_shifted{current_frame, current_bin} = ...
                horzcat(bin_fluo_data_shifted{current_frame, current_bin}, new_data.ms2{i}.Fluo(j));

        end;
    end;
    
    % Normalized slopes
    if ~isnan(normalized_slopes_array(i))
        bin_normalized_slopes{current_bin} = horzcat (bin_normalized_slopes{current_bin}, normalized_slopes_array(i));
    end;
    
end;



%% Calculating mean and STD

% Non-shifted data
bin_fluo_data_non_shifted_mean = zeros(size(bin_fluo_data_non_shifted));
bin_fluo_data_non_shifted_STD = zeros(size(bin_fluo_data_non_shifted));
bin_fluo_data_non_shifted_count = zeros(size(bin_fluo_data_non_shifted));
for i=1:size(bin_fluo_data_non_shifted,1)
    for j=1:size(bin_fluo_data_non_shifted,2)     
            bin_fluo_data_non_shifted_mean(i,j) = mean(bin_fluo_data_non_shifted{i,j});
            bin_fluo_data_non_shifted_STD(i,j) = std(bin_fluo_data_non_shifted{i,j});
            bin_fluo_data_non_shifted_count(i,j) = length(bin_fluo_data_non_shifted{i,j});
    end;
end;



% Shifted data
bin_fluo_data_shifted_mean = zeros(size(bin_fluo_data_shifted));
bin_fluo_data_shifted_STD = zeros(size(bin_fluo_data_shifted));
bin_fluo_data_shifted_count = zeros(size(bin_fluo_data_shifted));
bin_fluo_data_shifted_confidence_intervals = zeros([size(bin_fluo_data_shifted),2]);
bin_normalized_slopes_confidence_intervals = zeros(bins_count, 2);
bootfunc = @(x) mean(x);
for frame=1:1:size(bin_fluo_data_shifted,1)
    for bin=1:size(bin_fluo_data_shifted,2)
            bin_fluo_data_shifted_mean(frame,bin) = mean(bin_fluo_data_shifted{frame,bin});
            bin_fluo_data_shifted_STD(frame,bin) = std(bin_fluo_data_shifted{frame,bin});
            bin_fluo_data_shifted_count(frame,bin) = length(bin_fluo_data_shifted{frame,bin});
            
            
    end;
    
end;


% Normalized slopes
bin_normalized_slopes_mean = zeros(1, bins_count);
for bin = 1:bins_count 
    bin_normalized_slopes_mean(bin) = mean(bin_normalized_slopes{bin});
end;


display('Calculating bootstrap confidence intervals');

% Bootstrapping all bins

% Choosing the bin for boostrap plotting. Calculating bootstrap intervals
% for all bins
bootstrap_only_bin_number = max(find(bins_borders >= bootstrap_only_bin_value, 1, 'first')-1,1);


% Bootstrapping the data
for bin=1:size(bin_fluo_data_shifted,2)
    for frame=1:bootstrap_only_each_frame:size(bin_fluo_data_shifted,1)

            
            %% Calculating bootstrapped confidence intervals for the data
            if ... % bin==bootstrap_only_bin_number &&...   % Now bootstrapping all bins
               ~isempty(bin_fluo_data_shifted{frame,bin}) &&...
               length(bin_fluo_data_shifted{frame,bin})>1
                        bin_fluo_data_shifted_confidence_intervals(frame,bin,:) = bootci(bootstrap_samples_count, bootfunc, bin_fluo_data_shifted{frame,bin});
                        
            end;
            
            % Subtracting the mean if bootstrapping makes sense
            if bin_fluo_data_shifted_confidence_intervals(frame,bin,1)>0
                bin_fluo_data_shifted_confidence_intervals(frame,bin,1) = bin_fluo_data_shifted_mean(frame,bin) - bin_fluo_data_shifted_confidence_intervals(frame,bin,1);
            else
                bin_fluo_data_shifted_confidence_intervals(frame,bin,1) = 0;
            end;
            if bin_fluo_data_shifted_confidence_intervals(frame,bin,2)>0
                bin_fluo_data_shifted_confidence_intervals(frame,bin,2) = bin_fluo_data_shifted_confidence_intervals(frame,bin,2) - bin_fluo_data_shifted_mean(frame,bin);
            else
                bin_fluo_data_shifted_confidence_intervals(frame,bin,2) = 0;
            end;
            
            
    end;
    
    %% Bootstrapping the normalized slopes
    if ~isempty(bin_normalized_slopes{bin}) &&...
           length(bin_normalized_slopes{bin})>1
        bin_normalized_slopes_confidence_intervals(bin,:) = bootci(bootstrap_samples_count, bootfunc, bin_normalized_slopes{bin});

    end;
    
% % %     % Subtracting the mean if bootstrapping makes sense
% % %     if bin_normalized_slopes_confidence_intervals(bin,1)>0
% % %         bin_normalized_slopes_confidence_intervals(bin,1) = bin_normalized_slopes_mean(bin) - bin_normalized_slopes_confidence_intervals(bin,1);
% % %     else
% % %         bin_normalized_slopes_confidence_intervals(bin,1) = 0;
% % %     end;
% % %     if bin_normalized_slopes_confidence_intervals(bin,2)>0
% % %         bin_normalized_slopes_confidence_intervals(bin,2) = - bin_normalized_slopes_mean(bin) + bin_normalized_slopes_confidence_intervals(bin,2);
% % %     else
% % %         bin_normalized_slopes_confidence_intervals(bin,2) = 0;
% % %     end;
    
    
    fprintf('Bootstrapping progress: %.1f%%\n', bin/size(bin_fluo_data_shifted,2)*100);
end;


% % % % Subtracting the mean if bootstrapping makes sense
% % % bin_fluo_data_shifted_confidence_intervals(frame,bin,1) = bin_fluo_data_shifted_mean(frame,bin) - bin_fluo_data_shifted_confidence_intervals(frame,bin,1);
% % % bin_fluo_data_shifted_confidence_intervals(frame,bin,2) = bin_fluo_data_shifted_mean(frame,bin) - bin_fluo_data_shifted_confidence_intervals(frame,bin,2);

display('Bootstrap confidence intervals calculation completed');



%% Plotting the evolution of the mean fluorescence
figure(3);
clf;
subplot(1,2,1);
hold on;
subplot(1,2,2);
hold on;
str_legends=cell(1, bins_count);

for bin=1:bins_count
    % Non-shifted data
    subplot(1,2,1);
    plot(original_time_mesh,...
        bin_fluo_data_non_shifted_mean(:,bin)/fluo_per_polymerase, 'LineWidth', 2);
    
    % Shifted data
    subplot(1,2,2);
    if bin==bootstrap_only_bin_number
        
        error_handle = errorbar(time_mesh(1:bootstrap_only_each_frame:end),...
            bin_fluo_data_shifted_mean(1:bootstrap_only_each_frame:end,bin)/fluo_per_polymerase,...
            bin_fluo_data_shifted_confidence_intervals(1:bootstrap_only_each_frame:end, bin,1)/fluo_per_polymerase,...
            bin_fluo_data_shifted_confidence_intervals(1:bootstrap_only_each_frame:end, bin,2)/fluo_per_polymerase,...
        '-o', 'LineWidth', 2);
%         errorbar_tick(error_handle, 50);
    else
        plot(time_mesh,...
        bin_fluo_data_shifted_mean(:,bin)/fluo_per_polymerase, 'LineWidth', 2);
    end;
    
    % Saving the legend
    str_legends{bin} = sprintf('AP = %.2f', bins_borders(bin) + bin_width/2);
    
end;

% xlim_vector = [48, 60]; %[48, 56]
ylim_vector = [0,110];

subplot(1,2,1);
xlabel('Time, min');
ylabel('Number of polymerases');
xlim(xlim_vector);
ylim(ylim_vector);
title('Non-shifted data');
legend(str_legends, 'Location', 'northwest');

% Adding the fill
x_vertices = forced_start_nc_time + [0, 20, 20, 0];
y_vertices = [0,...
    k/(1+sqrt(l_min))^2 * (x_vertices(2) - x_vertices(1)),...
    k/(1+sqrt(l_max))^2 * (x_vertices(3) - x_vertices(1)),...
    0];


subplot(1,2,2);
xlabel('Time, min');
ylabel('Number of polymerases');
xlim(xlim_vector);
ylim(ylim_vector);
title('Shifted data');
legend(str_legends, 'Location', 'northwest');
fill(x_vertices, y_vertices, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'None');

% Adding the theoretical prediction to the second plot
subplot(1,2,2);
theor_time_mesh = forced_start_nc_time + (0:time_step:init_slope_length);
theor_slope_func = @(t) k/(1+sqrt(l))^2 * t; % k in min^{-1}, t in min
theor_slope_data = theor_slope_func(theor_time_mesh-forced_start_nc_time);

plot(theor_time_mesh, theor_slope_data, '--', 'LineWidth', 2, 'color', 'black');

% Adding the maximal particle number
temp_time_mesh = xlim_vector(1):time_step:xlim_vector(2);
plot(temp_time_mesh, ones(1,length(temp_time_mesh)) * L/(l+sqrt(l)), ':', 'color', 'black',...
    'LineWidth', 2);





% % % 
% % % %% Plotting the STD
% % % figure(4);
% % % clf;
% % % subplot(1,2,1);
% % % hold on;
% % % subplot(1,2,2);
% % % hold on;
% % % 
% % % for bin=1:bins_count
% % %     % Non-shifted data
% % %     subplot(1,2,1);
% % %     plot(original_time_mesh,...
% % %         bin_fluo_data_non_shifted_STD(:,bin)/fluo_per_polymerase, 'LineWidth', 2);
% % %     
% % %     % Shifted data
% % %     subplot(1,2,2);
% % %     plot(new_time_mesh,...
% % %         bin_fluo_data_shifted_STD(:,bin)/fluo_per_polymerase, 'LineWidth', 2);
% % % end;
% % % 
% % % % xlim_vector = [48, 60]; % 56
% % % ylim_vector = [0,32];
% % % 
% % % subplot(1,2,1);
% % % xlabel('Time, min');
% % % ylabel('STD of N');
% % % xlim(xlim_vector);
% % % ylim(ylim_vector);
% % % title('Non-shifted data');
% % % 
% % % subplot(1,2,2);
% % % xlabel('Time, min');
% % % ylabel('STD of N');
% % % xlim(xlim_vector);
% % % ylim(ylim_vector);
% % % title('Shifted data');



%% Plotting variance as a function of mean as suggested by Tkacik
figure(7);
clf;

for bin=1:bins_count
    % Non-shifted data
    subplot(1,2,1);
    scatter(bin_fluo_data_non_shifted_mean(:,bin)/fluo_per_polymerase,...
        (bin_fluo_data_non_shifted_STD(:,bin)/fluo_per_polymerase).^2, 'LineWidth', 2);
    hold on;
    
    % Shifted data
    subplot(1,2,2);
    scatter(bin_fluo_data_shifted_mean(:,bin)/fluo_per_polymerase,...
        (bin_fluo_data_shifted_STD(:,bin)/fluo_per_polymerase).^2, 'LineWidth', 2);
    hold on;
end;

ylim_vector = [0, 500];

subplot(1,2,1);
xlabel('<N>');
ylabel('Var(N)');
% xlim(xlim_vector);
ylim(ylim_vector);
title('Non-shifted data');
hold off;

subplot(1,2,2);
xlabel('<N>');
ylabel('Var(N)');
% xlim(xlim_vector);
ylim(ylim_vector);
title('Shifted data');
hold off;



%% Plotting the number of spots in a bin
figure(5);
clf;
subplot(1,2,1);
hold on;
subplot(1,2,2);
hold on;

for bin=1:bins_count
    % Non-shifted data
    subplot(1,2,1);
    plot(original_time_mesh,...
        bin_fluo_data_non_shifted_count(:,bin), 'LineWidth', 2);
    
    % Shifted data
    subplot(1,2,2);
    plot(time_mesh,...
        bin_fluo_data_shifted_count(:,bin), 'LineWidth', 2);
end;

% xlim_vector = [48, 60]; % 56
% % ylim_vector = [0,25];

subplot(1,2,1);
xlabel('Time, min');
ylabel('Number of data points');
xlim(xlim_vector);
% ylim(ylim_vector);
title('Non-shifted data');

subplot(1,2,2);
xlabel('Time, min');
ylabel('Number of data points');
xlim(xlim_vector);
% ylim(ylim_vector);
title('Shifted data');





% 
% %% Looking how these fluctuations are distributed across the ensemble
% figure(6);
% clf;
% % % subplot(1,2,1);
% % % hold on;
% % % subplot(1,2,2);
% % % hold on;
% 
% selected_new_mesh_frame_index = 847;
% % hold on;
% for bin=1:bins_count
%     subplot(1, bins_count, bin);
%     histogram(bin_fluo_data_shifted{selected_new_mesh_frame_index, bin}/fluo_per_polymerase, 1:8:80,...
%         'Normalization', 'count');
%     title(sprintf('AP = %.2f. T = %.0f min', ...
%         bins_borders(bin) + bin_width/2,...
%         new_time_mesh(selected_new_mesh_frame_index)));
%     xlabel('N');
%     ylabel('PDF');
% end;



%% Plotting data distribution in bins
% figure(8);




%% Now I want to check how changes with time the number of data I calculate
% means and stds over

%% Outputtingthe AP interval I consider
fprintf('AP interval for %s: [%.2f, %.2f]\n', dataset_name, min_ms2_AP, max_ms2_AP);


%% Let's find the maximal slope for the shifted data. 
% For this let's first find all slopes. 

% plot(new_time_mesh,...
%         bin_fluo_data_shifted_mean(:,bin)/fluo_per_polymerase, 'LineWidth', 2);

bin_mean_slopes = zeros(1, bins_count);

% Selecting the right time interval
time_interval_slope = 2;    % where to look for the slope, in mins
start_time_mins = forced_start_nc_time;
t_ind_slope_start = find(time_mesh > start_time_mins, 1, 'first');
t_ind_slope_end = find(time_mesh > start_time_mins+time_interval_slope, 1, 'first');

indices_to_consider = t_ind_slope_start:t_ind_slope_end;

for bin = 1:bins_count
    
    % From the indices to consider I have to choose only those that are not
    % NaN
    indices_not_nan = find(~isnan(bin_fluo_data_shifted_mean(:,bin)));
        
    indices_to_consider_not_nan = intersect(indices_to_consider, indices_not_nan);
    
    
    fit_result= ((time_mesh(indices_to_consider_not_nan)' - start_time_mins)\...
        [ones(length(indices_to_consider_not_nan),1), bin_fluo_data_shifted_mean(indices_to_consider_not_nan,bin)]);
    bin_mean_slopes(bin) = fit_result(2);
end;

bin_mean_slopes = bin_mean_slopes/fluo_per_polymerase;

% Now print out the maximal slope over all bins
max_experimental_slope = max(bin_mean_slopes);
theoretical_slope = k/(1+sqrt(l))^2;
fprintf('\nMaximal experimental slope: %.2f pol/min\n', max_experimental_slope);
fprintf('Theoretical slope: %.2f pol/min\n', theoretical_slope);
fprintf('Percent of the theoretical slope: <=%.0f%%\n', max_experimental_slope/theoretical_slope*100);


%% Outputting the processed data
output_folder = './processed_data/';
output_filename = sprintf('%s_%s_nc_%i.mat', gene_name, dataset_name, nuc_cyc);
output_full_path = strcat(output_folder, output_filename);

save(output_full_path, 'bin_fluo_data_shifted_mean', 'new_time_mesh', 'bins_borders','bins_count', 'bins_centers', 'bin_width', ...
        'bin_fluo_data_shifted_confidence_intervals', 'bin_fluo_data_shifted_STD', 'bin_fluo_data_shifted', 'bin_fluo_data_shifted_count',...
        'slopes_array', 'normalized_slopes_array', 'bin_normalized_slopes', ...
        'bin_normalized_slopes_mean', 'bin_normalized_slopes_confidence_intervals', 'ms2_combined',...
        'ms2_spots_count', 'forced_start_nc_time', 'intersct_array');



%% Old plot functions
% plot_fig_1b;
% plot_fig_3;
% plot_fig_4;   
% plot_slopes_vs_AP;
% plot_poll_number_vs_AP;
    

    end;
end;

%% Plotting figures for the article
plot_mean_evolution_article;
plot_slopes_histogram_article;
plot_slopes_vs_AP_article_1_plot;








