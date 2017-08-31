
set (0, 'DefaultAxesFontSize', 20);
clear;



%% Constants
constants;



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
for i=1:ms2_count
    plot (ms2_combined(i).Frame * mins_per_frame, ms2_combined(i).Fluo/fluo_per_polymerase);
end;

% Label plot
xlabel('Time, min');
ylabel('Number of active polymerases');
title('Roughly adjusted fluorescence data');
hold off;
    


%% Group different spots into AP bins

% Get the average AP position of each trace
ms2_mean_AP = zeros(1,ms2_count);
for i=1:ms2_count
    ms2_mean_AP(i) = mean(ms2_combined(i).APpos);
end;

% If the AP limits are not manually defined, use the min and max of the observed values
if ~manual_bins
    min_ms2_AP = min(ms2_mean_AP);
    max_ms2_AP = max(ms2_mean_AP);
end;

% Define bin locations
bins_borders = (min_ms2_AP-bin_width/2):bin_width:(max_ms2_AP + bin_width*3/2);
bins_count = length(bins_borders)-1;
bins_centers = (bins_borders(1:end-1) + bins_borders(2:end))/2;

% Assign each trace to a bin
ms2_index_to_bin_table = zeros(1, ms2_count);
for i=1:ms2_count
    current_bin = ceil((ms2_mean_AP(i)-bins_borders(1))/bin_width);
    ms2_index_to_bin_table (i) = current_bin;
end;



%% Calculate max frame number (maximum time) for the collected traces
max_frame_number = max(ms2_combined(1).Frame);
for i = 2:ms2_count
    max_frame_number = max(max_frame_number, max(ms2_combined(i).Frame));
end;
min_frame_number = 1;



%% Let's now try to determine initial slope of each of the curves so as to be able to sync them
% New feature: now the initial slope will not be identified on a fixed
% interval, but rather from the beginning of the nc to the max.
% luminescence value achieved. So slope_end_frame is floating

% slope_start_frame = round(init_slope_interval(1)/mins_per_frame);
% slope_end_frame = round(init_slope_interval(2)/mins_per_frame);



% For each ms2 spots extract the data corresponding to this time, and fit
% it with a straight line
intersct_array = zeros(1, ms2_count);
slopes_array = zeros(1, ms2_count);
coefs_array = zeros(ms2_count, 2);
for i =1: ms2_count
    points_count = length(ms2_combined(i).Frame);
    
    %% Detecting the first frame of the nc
    current_fluo_cut = ms2_combined(i).Fluo;
    % Removing all contributions which lie before the current nc
    current_fluo_cut (ms2_combined(i).Frame < forced_start_nc_time/mins_per_frame) = 0;
    slope_start_index = find(current_fluo_cut > 0, 1);
    if isempty(slope_start_index)
        continue;
    end;
    
    %% Finding the frame where the slope ends.
    % I am working under the assumption that maxima from other nc do not
    % interfere (that ncs are correctly labeled)
% %     slope_end_index = find(ms2_combined(i).Fluo == max(ms2_combined(i).Fluo));

    slope_end_index = slope_start_index + ceil(init_slope_length / mins_per_frame);
    slope_end_index = min (slope_end_index, points_count);
%     disp(slope_start_index);
%     disp(slope_end_index);

    %% Selecting points corresponding to the slope
    selected_fluo_points = [];
    selected_time_points = [];
    selected_points = [];
    count = 0;
    for j=1:points_count
        if j >= slope_start_index && j <= slope_end_index
            count = count+1;
            selected_points (count) = j;
            selected_time_points (count) = ms2_combined(i).Frame(j) * mins_per_frame;
            selected_fluo_points (count) = ms2_combined(i).Fluo(j);
        end;
    end;
    
    % Fitting with a straight line
    if count > 1
        coefs = polyfit(selected_time_points, selected_fluo_points, 1);
        coefs_array (i, :) = coefs;
        % Calculating the time of the beginning of the n. cycle
        intersct = - coefs(2)/coefs(1);
        intersct_array(i) = intersct;
        slopes_array (i) = coefs(1)/fluo_per_polymerase;    % These slopes will be in pols/min
    end;
end;









% Filtering out the curves which in the interval of interest are below a
% threshold value
%         time_indi
%         fluo_in_the_analyzed_interval = (ms2_combined(ms2_count).Fluo(:)
%         if max(ms2_combined(ms2_count).Fluo(:)/fluo_per_polymerase) < N_filter_threshold
%             continue;
%         end;





%% Now plotting the data again, but shifting the data for each curve in such a way
% that the beginning of nc14 is at t=50mins

% Plotting all together for nc14 only


figure(2);
clf;
hold on;
% xlim(xlim_vector);
for i=1:ms2_count
    if ms2_combined(i).nuc_cyc ~= nuc_cyc; continue; end;
    
    % Taking only those ms2 for which the calculated intersect was > 0
    if intersct_array(i)>0
        plot ((ms2_combined(i).Frame * mins_per_frame) + (forced_start_nc_time - intersct_array(i)), ...
            ms2_combined(i).Fluo/fluo_per_polymerase);
        hold on;
        
%         %% Plotting the slope
%         x_data_mins = forced_start_nc_time + (0:0.05:12);
%         y_data_pols = (coefs_array(i, 1) .* (x_data_mins - forced_start_nc_time)) / fluo_per_polymerase;
%         plot(x_data_mins, y_data_pols, 'k--');
        
        xlim(xlim_vector);
        ylim([0,100]);
%         hold off;
    end;

end;
xlabel('Time, min');
ylabel('Number of active polymerases');

% xlim([47, 70]);
% hold off;
title('Shifted fluorescence data');

% Adding the theoretical prediction
% hold on;
theor_time_mesh = forced_start_nc_time:0.05:(forced_start_nc_time+4);
theor_slope_func = @(t) k/(1+sqrt(l))^2 * t; % k in min^{-1}, t in min
theor_slope_data = theor_slope_func(theor_time_mesh-forced_start_nc_time);
plot(theor_time_mesh, theor_slope_data, '--', 'LineWidth', 2, 'color', 'black');


% Adding the maximal particle number
temp_time_mesh = xlim_vector(1):0.05:xlim_vector(2);
plot(temp_time_mesh, ones(1,length(temp_time_mesh)) * L/(l+sqrt(l)), ':', 'color', 'black',...
    'LineWidth', 2);
hold off;


%%% ===
if bl_stop_early
    return;
end;
%%% ==


%% Creating a new, more refined, time mesh; shiftind the data and projecting the shifted data on 
% the new time mesh

% Creating a new, more refined, time mesh for the shifted data
refining_factor = 10;
original_time_mesh = (1:max_frame_number)*mins_per_frame;
new_time_step = mins_per_frame/refining_factor;    % in mins
new_time_mesh = original_time_mesh(1): new_time_step :original_time_mesh(end);
new_time_mesh_length = length(new_time_mesh);

% Initializing new data array
new_data = {};
new_data.ms2 = cell(1, ms2_count);

% Parsing data points
for i=1:ms2_count
    % Checking if we have a meaningful time shift, and skipping if not
    current_shift = forced_start_nc_time - intersct_array(i);
    % Conditions: 1. Intersect time should be positive
    % 2. The time shift should be withing the time_shift_limit (10s)
    % 3. The slope should be positive
    if ~(intersct_array(i) > 0 && abs(current_shift) <= time_shift_limit && slopes_array(i)>0)
        % Also, removing the slopes that don't fit this criteria
        slopes_array(i) = NaN;
        continue;
    end;
    
    original_time_points = ms2_combined(i).Frame * mins_per_frame;  % in mins
    original_fluo = ms2_combined(i).Fluo;
    
    shifted_time_points = original_time_points + current_shift;
    
   
    % Determining new frame numbers on which to project the data points
    query_time_start = min(new_time_mesh(new_time_mesh>=shifted_time_points(1)));
    query_new_frames_ind_start = find(new_time_mesh == query_time_start);
    query_time_end = max(new_time_mesh(new_time_mesh<=shifted_time_points(end)));
    query_new_frames_ind_end = find(new_time_mesh == query_time_end);
    query_time_points = query_time_start: new_time_step :query_time_end;
    
    % Getting interpolated fluorescence data on the new time mesh frames
    fluo_at_query_time_points = interp1(shifted_time_points, original_fluo,...
        query_time_points, 'linear');
    
    % Saving results to new data cell array
    new_data.ms2{i}.newFrame = query_new_frames_ind_start:query_new_frames_ind_end;
    new_data.ms2{i}.Time = query_time_points;
    new_data.ms2{i}.Fluo = fluo_at_query_time_points;
    
end;


%% Normalizing the slope to the theoretical maximum
theoretical_slope = k/(1+sqrt(l))^2;
normalized_slopes_array = slopes_array / theoretical_slope;




%% Now plotting non-shifted and shifted data grouped in bins

% Initializing
% Non-shifted data
% frame_bins_count = ceil(max_frame_number-min_frame_number+1);

bin_fluo_data_non_shifted = cell(max_frame_number, bins_count);
% bin_fluo_data_count_non_shifted = cell(max_frame_number, bins_count);
% Shifted data
bin_fluo_data_shifted = cell(new_time_mesh_length, bins_count);
bin_normalized_slopes = cell(1, bins_count);
% bin_fluo_data_count_shifted = cell(new_time_mesh_length, bins_count);


% Collecting data in bins
for i = 1: ms2_count
    current_bin = ms2_index_to_bin_table(i);
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
        
        error_handle = errorbar(new_time_mesh(1:bootstrap_only_each_frame:end),...
            bin_fluo_data_shifted_mean(1:bootstrap_only_each_frame:end,bin)/fluo_per_polymerase,...
            bin_fluo_data_shifted_confidence_intervals(1:bootstrap_only_each_frame:end, bin,1)/fluo_per_polymerase,...
            bin_fluo_data_shifted_confidence_intervals(1:bootstrap_only_each_frame:end, bin,2)/fluo_per_polymerase,...
        '-o', 'LineWidth', 2);
%         errorbar_tick(error_handle, 50);
    else
        plot(new_time_mesh,...
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
theor_time_mesh = forced_start_nc_time + (0:new_time_step:init_slope_length);
theor_slope_func = @(t) k/(1+sqrt(l))^2 * t; % k in min^{-1}, t in min
theor_slope_data = theor_slope_func(theor_time_mesh-forced_start_nc_time);

plot(theor_time_mesh, theor_slope_data, '--', 'LineWidth', 2, 'color', 'black');

% Adding the maximal particle number
temp_time_mesh = xlim_vector(1):new_time_step:xlim_vector(2);
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
    plot(new_time_mesh,...
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
t_ind_slope_start = find(new_time_mesh > start_time_mins, 1, 'first');
t_ind_slope_end = find(new_time_mesh > start_time_mins+time_interval_slope, 1, 'first');

indices_to_consider = t_ind_slope_start:t_ind_slope_end;

for bin = 1:bins_count
    
    % From the indices to consider I have to choose only those that are not
    % NaN
    indices_not_nan = find(~isnan(bin_fluo_data_shifted_mean(:,bin)));
        
    indices_to_consider_not_nan = intersect(indices_to_consider, indices_not_nan);
    
    
    fit_result= ((new_time_mesh(indices_to_consider_not_nan)' - start_time_mins)\...
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








