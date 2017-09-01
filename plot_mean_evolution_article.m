%% Plotting the mean fluorescence evolution curve for each gene



%% Constants
constants;

% Local constants
fill_opacity = 0.2;
x_lim_vec = [0, 15];
y_lim_vec = [0, 60];
AP_values_for_plot_array = [0.2, 0.3, 0.45];
take_every_which_curve = 1;
curves_per_construct = 100;
% fig_side = 350;
% fig_side_ratio = [1.25, 1];

% bl_save_figures = true;
% bl_save_figures = false;

local_color_sequence = [
	0    0.4470    0.7410;... % blue
    0.9290    0.6940    0.1250;... % orange/yellow
    0.8500    0.3250    0.0980... % red
	];
local_curve_style_sequence = {'-', '--', ':'};

% slope_curve_color = [139  69  19] / 255;	% Saddle brown
slope_curve_color = [0 128 128] / 255;	% Teal



%% Initialize
AP_values_number = length(AP_values_for_plot_array);
max_curve_number = length(dataset_names_array) * length(gene_names_array) * 2;
saved_curves_colors = zeros(max_curve_number, 3);

% Initialize figure
fig_hand_nc13 = figure(4);
fig_size = fig_side * fig_side_ratio;
set_article_figure_size (fig_hand_nc13, 1, 1, 1);
clf;
hold on;



%% To be able to calculate the average curve, we have to collect all meshes
% I am using the same mesh for all data sets, so no need to store different meshes.




%% Plot
total_plotted_curves_counter = 0;
fluo_data_combined = [];
% time_mesh_length = length(time_mesh);	% Remove this line
% average_evolution_profile = zeros(time_mesh_length, 1);
figure(fig_hand_nc13);
hold on;
for nuc_cyc = [13, 14]
    if nuc_cyc == 13
        time_shift = forced_start_nc_time_13;
    else
        time_shift = forced_start_nc_time_14;
    end;

    %% Parsing different genes first to plot the mean values
    for gene_index = 1:length(gene_names_array)
        gene_name = gene_names_array{gene_index};
%         construct_curves_plotted = 0;
        for construct_index = 1:length(dataset_names_array)
            dataset_name = dataset_names_array{construct_index};
            % Skipping Knirps with no_shadow nc 13, because for it we have no data
            if nuc_cyc ==13 && gene_index == 3 && construct_index == 3
                continue;
            end;
			
            % Load data
            input_filename = sprintf('%s_%s_nc_%i.mat', gene_name, dataset_name, nuc_cyc);
            input_full_path = strcat(input_folder, input_filename);
            load (input_full_path);
			
            % Average over bins, keep frames
            fluo_bin_average = mean(bin_fluo_data_mean, 1, 'omitnan') / fluo_per_polymerase;
 
			h_chart = plot(time_mesh - time_shift, fluo_bin_average, 'color', local_color_sequence(gene_index, :), 'LineWidth', 1);
			
			% Store data for the total correctly-weighted average.
			% I drop here the bin index and reshape the array.
			% Permute indices setting bins at the end: [bin, frame, trace] -> [frame, trace, bin]
			tmp_data = permute(bin_fluo_data_new, [2, 3, 1]);
			% Drop the last index to get rid of the bin dimension -> [frame, trace]
			tmp_data = reshape(tmp_data, time_mesh_length, []);
			
			
			
			%% Need to shift data, so that the start of the nuclear cycle be at 0
			% Determine shift in frames
			time_shift_frames = round((time_shift - new_time_start) / time_step);	% Check later if new_time_start is placed correctly
			% Need to cut off first frames of the array to shift.
			tmp_data = tmp_data(time_shift_frames + 1 : end, :);
			% Restore the length with NaN
			nan_array = ones(time_shift_frames, size(tmp_data, 2)) * NaN;
			tmp_data = cat(1, tmp_data, nan_array);
			
			% Concatenate with earlier stored data along the trace dimension (2)
			fluo_data_combined = cat(2, fluo_data_combined, tmp_data);
			
% %                     average_evolution_profile = average_evolution_profile + cur_data;
% 			total_plotted_curves_counter = total_plotted_curves_counter + 1;
% 			saved_meshes_and_curves_and_weigths{total_plotted_curves_counter} = ...
% 				[time_mesh - time_shift; fluo_bin_average'; sum(bin_fluo_data_shifted_count, 2)'];
%                 saved_curve_weight(total_plotted_curves_counter, 1) = sum(bin_fluo_data_shifted_count, 2);  % We sum the weight over the bins
%             end;
        end;
    end;

end;



%% Add theoretical slope nc 13
time_step = 5/60;		% Remove
figure(fig_hand_nc13);
% Maximal slope
forced_start_nc_time = 0;
theor_line_width = 2;
theor_time_mesh = forced_start_nc_time : time_step : 2*init_slope_length;
theor_slope_func = @(t) k/(1+sqrt(l))^2 * t; % k in min^{-1}, t in min
theor_slope_data = theor_slope_func(theor_time_mesh - forced_start_nc_time);
hnd_theor_max = plot(theor_time_mesh, theor_slope_data, '--', 'LineWidth', theor_line_width, 'color', slope_curve_color);
% Slow slope
theor_slope_func = @(t) (k/60*tau_exp - 1)/ (tau_exp/60) / (k/60 * tau_exp + l -1) * t; % k in min^{-1}, t in min
theor_slope_data = theor_slope_func(theor_time_mesh-forced_start_nc_time);
hnd_theor_abortive = plot(theor_time_mesh, theor_slope_data, '-', 'LineWidth', theor_line_width, 'color', slope_curve_color);
% Max number of the original l-TASEP.
x_lim = xlim();
theor_max_number = L/(l+sqrt(l));
plot(x_lim, theor_max_number * ones(1,2), '--', 'LineWidth', theor_line_width, 'color', 'black');


% % % %% Adding a theoretical maximum prediction of the modified model.
% % % theor_max_with_pause = L / (k/60 * tau_exp + l - 1);
% % % x_lim_vec = xlim();
% % % theor_data = ones(1, 2) .* theor_max_with_pause;
% % % plot(x_lim_vec, theor_data, 'k', 'LineWidth', theor_line_width);


%% Plot the average trace
mean_trace = mean(fluo_data_combined, 2, 'omitnan');
% Convert to polymerases
mean_trace = mean_trace / fluo_per_polymerase;

% % % 
% % % t_start = 0;
% % % t_end = 15; % in mins
% % % local_t_mesh = t_start:new_time_step:t_end;
% % % local_t_mesh_length = length(local_t_mesh);
% % % % Interpolating and getting data from all curves at these time steps
% % % organized_mean_curves_data = zeros(total_plotted_curves_counter, local_t_mesh_length);
% % % weights = zeros(total_plotted_curves_counter, local_t_mesh_length);
% % % for i = 1:total_plotted_curves_counter
% % %     weights(i, :) = interp1(saved_meshes_and_curves_and_weigths{i}(1, :),...
% % %         saved_meshes_and_curves_and_weigths{i}(3, :), local_t_mesh);
% % %     organized_mean_curves_data(i, :) = interp1(saved_meshes_and_curves_and_weigths{i}(1, :),...
% % %         saved_meshes_and_curves_and_weigths{i}(2, :), local_t_mesh);
% % % end;
% % % % I also have to take into accounts weights of individual already averaged
% % % % curves
% % % 
% % % % Average
% % % nan_indices = find(isnan(organized_mean_curves_data));
% % % weights(nan_indices) = 0;
% % % organized_mean_curves_data(nan_indices) = 0;
% % % 
% % % 
% % % mean_mean_curve = sum(organized_mean_curves_data .* weights, 1) ...
% % %     ./ sum(weights, 1);
% Plot
plot(time_mesh, mean_trace, '-', 'color', 'k', 'LineWidth', 3);



%% Reorder curves
uistack(hnd_theor_max, 'bottom');
uistack(hnd_theor_abortive, 'bottom');


%% Adjusting plot parameters
% NC 13
figure(fig_hand_nc13);
xlim(x_lim_vec);
ylim(y_lim_vec);
xlabel('Time since the start of a nuc. cyc. (min)'); 
ylabel('N');
grid on;
box on;

%% Saving the figure
if bl_save_figures
    % Adjusting size
    fig_hand_nc13.PaperPositionMode = 'auto';
    fig_hand_nc13.Units = 'Inches';
    pos = fig_hand_nc13.Position;
    set(fig_hand_nc13, 'PaperUnits','Inches','PaperSize', [pos(3), pos(4)]);
% % %     % PNG
% % %     output_filename_svg = sprintf('02_Mean_evolution_article.png');
% % %     output_full_path_svg = strcat(output_figures_folder, output_filename_svg);
% % %     display('Saving figure');
% % %     fig_hand_nc13.PaperPositionMode = 'auto';
% % %     saveas(fig_hand_nc13, output_full_path_svg, 'png');
% % %     display('Figure saved');
    % PDF
    output_filename_pdf = sprintf('02_Mean_evolution_article.pdf');
    output_full_path_pdf = strcat(output_figures_folder, output_filename_pdf);
    display('Saving figure');
    print(fig_hand_nc13, output_full_path_pdf, '-dpdf', '-r0');
    display('Figure saved');
end;


