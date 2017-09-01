



%% Plotting the maximal fluorescence curve for each gene





%% Constants
dataset_names_array = {'bac', 'no_primary', 'no_shadow'};
dataset_short_names_array = {'bac', 'no_pr', 'no_sh'};
gene_names_array = {'HunchBack', 'Knirps', 'SNAIL'};
gene_short_names_array = {'Hb', 'Kn', 'Sn'};
input_folder = './processed_data/';
fill_opacity = 0.2;
x_lim_vector = [34, 65];
y_lim_vector = [0, 110];
AP_values_for_plot_array = [0.2, 0.3, 0.45];
hist_color = [117, 145, 41] / 255;
% image_width = 500; % in px (?)
% image_sides_ratio = 0.7;
% abortive_injection_slope_modifier = 0.34;
take_every_which_curve = 1;
curves_per_construct = 100;
% bl_save_figures = true;
% bl_save_figures = false;

%% Initialization
AP_values_number = length(AP_values_for_plot_array);
max_curve_number = length(dataset_names_array) * length(gene_names_array) * 2;
saved_curves_colors = zeros(max_curve_number, 3);
local_color_sequence = [0    0.4470    0.7410;
                        0.9290    0.6940    0.1250;
                        0.8500    0.3250    0.0980];
local_curve_style_sequence = {'-', '--', ':'};



% Initializing the figure
fig_hand_nc13 = figure(4);
fig_size = 350 * [8/7, 1];
set_my_fig_size (fig_hand_nc13, fig_size);
clf;
hold on;
% fig_hand_nc14 = figure(5);
% set_my_fig_size (fig_hand_nc14, image_width * [1, image_sides_ratio]);
% clf;
% hold on;

%% To be able to calculate the average behavior, we have to collect all meshes
% and all data since meshes are different for each mean curve
saved_meshes_and_curves_and_weigths = {};
% saved_curve_weight_by_frame = {};


%% Plot
total_plotted_curves_counter = 0;
average_evolution_profile = zeros(new_time_mesh_length, 1);
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
            % Skipping Knirps with no_shadow nc 13, because it we have no data
            if nuc_cyc ==13 && gene_index == 2 && construct_index == 3
                continue;
            end;
            % Loading data
            input_filename = sprintf('%s_%s_nc_%i.mat', gene_name, dataset_name, nuc_cyc);
            input_full_path = strcat(input_folder, input_filename);
            load (input_full_path);
            % Averaging over bins
            fluo_bin_average = mean(bin_fluo_data_shifted_mean(:, :), 2, 'omitnan') / fluo_per_polymerase;
 
% % %             for bin_ind = 1:9 %bins_count
% % %             %     subplot(AP_values_number, 1, AP_ind);
% % %             %     hold on;
% % %         %         AP_value_for_plot = AP_values_for_plot_array(bin_ind);
% % %             %% Initialization
% % %             bin_index_for_plot = bin_ind; % ceil((AP_value_for_plot-bins_borders(1))/bin_width);
% % %             bl_curve_empty = zeros(1, max_curve_number);
% % %             % Counting not NaN indices
% % %             number_of_data_points = sum(~isnan(bin_fluo_data_shifted_mean(:, bin_index_for_plot)));
% % % % %             % Assembling the curve style
% % % % %             curve_style = local_curve_style_sequence{construct_index};
% % %             % Averaging over AP
% % %             % Plotting the bin we found
% % %             if number_of_data_points
% % %                 % Do not plot each curve
% % % %                     bl_plot = mod(counter, take_every_which_curve) == 0;
% % %                   bl_plot = construct_curves_plotted < curves_per_construct;
% % %                 
% % %                 
% % %             end;
            % Plotting the averaged bin data
%             if bl_plot
% % %                     cur_data = bin_fluo_data_shifted_mean(:, bin_index_for_plot)/fluo_per_polymerase;
                h_chart = plot(new_time_mesh - time_shift, fluo_bin_average,...
                     'color', local_color_sequence(gene_index, :),...
                    'LineWidth', 1);
%                 construct_curves_plotted = construct_curves_plotted + 1;
                % Preparing the average data
%                     average_evolution_profile = average_evolution_profile + cur_data;
                total_plotted_curves_counter = total_plotted_curves_counter + 1;
                saved_meshes_and_curves_and_weigths{total_plotted_curves_counter} = ...
                    [new_time_mesh - time_shift; fluo_bin_average'; sum(bin_fluo_data_shifted_count, 2)'];
%                 saved_curve_weight(total_plotted_curves_counter, 1) = sum(bin_fluo_data_shifted_count, 2);  % We sum the weight over the bins
%             end;
        end;
    end;

end;

%% Adding theoretical slope nc 13
figure(fig_hand_nc13);
% Maximal slope
forced_start_nc_time = 0;
theor_line_width = 2;
theor_time_mesh = forced_start_nc_time:new_time_step:x_lim_vector(2);
theor_slope_func = @(t) k/(1+sqrt(l))^2 * t; % k in min^{-1}, t in min
theor_slope_data = theor_slope_func(theor_time_mesh-forced_start_nc_time);
plot(theor_time_mesh, theor_slope_data, '--', 'LineWidth', theor_line_width, 'color', 'black');
% Slow slope
theor_slope_func = @(t) (k/60*tau_exp - 1)/ (tau_exp/60) / (k/60 * tau_exp + l -1) * t; % k in min^{-1}, t in min
theor_slope_data = theor_slope_func(theor_time_mesh-forced_start_nc_time);
plot(theor_time_mesh, theor_slope_data, '-', 'LineWidth', theor_line_width, 'color', 'black');
% Max number
x_lim = xlim();
theor_max_number = L/(l+sqrt(l));
plot(x_lim, theor_max_number * ones(1,2), '--', 'LineWidth', theor_line_width, 'color', 'black');


%% Adding a theoretical maximum prediction
theor_max_with_pause = L / (k/60 * tau_exp + l - 1);
x_lim_vec = xlim();
theor_data = ones(1, 2) .* theor_max_with_pause;
plot(x_lim_vec, theor_data, 'k', 'LineWidth', theor_line_width);


%% Calculating the weigthed mean curve
t_start = 0;
t_end = 15; % in mins
local_t_mesh = t_start:new_time_step:t_end;
local_t_mesh_length = length(local_t_mesh);
% Interpolating and getting data from all curves at these time steps
organized_mean_curves_data = zeros(total_plotted_curves_counter, local_t_mesh_length);
weights = zeros(total_plotted_curves_counter, local_t_mesh_length);
for i = 1:total_plotted_curves_counter
    weights(i, :) = interp1(saved_meshes_and_curves_and_weigths{i}(1, :),...
        saved_meshes_and_curves_and_weigths{i}(3, :), local_t_mesh);
    organized_mean_curves_data(i, :) = interp1(saved_meshes_and_curves_and_weigths{i}(1, :),...
        saved_meshes_and_curves_and_weigths{i}(2, :), local_t_mesh);
end;
% I also have to take into accounts weights of individual already averaged
% curves

% Average
nan_indices = find(isnan(organized_mean_curves_data));
weights(nan_indices) = 0;
organized_mean_curves_data(nan_indices) = 0;


mean_mean_curve = sum(organized_mean_curves_data .* weights, 1) ...
    ./ sum(weights, 1);
% Plot
plot(local_t_mesh, mean_mean_curve, 'color', 'k', 'LineWidth', 3);



%% Adjusting plot parameters
% NC 13
figure(fig_hand_nc13);
x_lim_vec = [0, 15];
xlim(x_lim_vec);
ylim([0, 60]);
xlabel('Time since the start of a nuc. cyc. (min)'); 
ylabel('N');
grid on;
box on;

%% Saving the figure
if bl_save_figures
    % Adjusting size
    fig_hand_nc13.PaperPositionMode = 'auto';
    set(fig_hand_nc13,'Units','Inches');
    pos = get(fig_hand_nc13,'Position');
    set(fig_hand_nc13,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    % PNG
    output_filename_svg = sprintf('02_Mean_evolution_article.png');
    output_full_path_svg = strcat(output_figures_folder, output_filename_svg);
    display('Saving figure');
    fig_hand_nc13.PaperPositionMode = 'auto';
    saveas(fig_hand_nc13, output_full_path_svg, 'png');
    display('Figure saved');
    % PDF
    output_filename_pdf = sprintf('02_Mean_evolution_article.pdf');
    output_full_path_pdf = strcat(output_figures_folder, output_filename_pdf);
    display('Saving figure');
    print(fig_hand_nc13, output_full_path_pdf, '-dpdf', '-r0');
    display('Figure saved');
end;

% % %     %% Adding the legend
% % %     legend_size = sum(~bl_curve_empty);
% % %     legend_string = cell(legend_size, 1);
% % %     curves_counter = 1;
% % %     legend_counter = 1;
% % %     for nuc_cyc = [13, 14]
% % %         for gene_index = 1:length(gene_names_array)
% % %             for i = 1:length(dataset_names_array)
% % %                 if ~bl_curve_empty(curves_counter)
% % %                     legend_string{legend_counter} = sprintf('nc %i, %s, %s', nuc_cyc, gene_short_names_array{gene_index}, dataset_short_names_array{i});
% % %                     legend_string{legend_counter} = strrep(legend_string{legend_counter}, '_', ' ');
% % %                     legend_counter = legend_counter + 1;
% % %                 end;
% % %                 curves_counter = curves_counter + 1;
% % %             end;
% % %         end;
% % %     end;
% % % %     if AP_ind == 1
% % % %         legend(legend_string, 'Location', 'northeast', 'FontSize', 12);
% % % %     end;


    %% Now parsing again to plot the fill (for the same AP)
    % Two for loops are required for correct legends
%     for AP_ind = 1:AP_values_number
%     subplot(AP_values_number, 1, AP_ind);
% %     hold on;
% 
%     AP_value_for_plot = AP_values_for_plot_array(AP_ind);



% % %     %% Initialization
% % % %     bin_index_for_plot = ceil((AP_value_for_plot-bins_borders(1))/bin_width);
% % %     counter = 1;
% % %     for gene_index = 1:length(gene_names_array)
% % %         gene_name = gene_names_array{gene_index};
% % %         for construct_index = 1:length(dataset_names_array)
% % %             dataset_name = dataset_names_array{construct_index};
% % %         
% % %             for nuc_cyc = [13, 14]
% % % 
% % %                 if ~bl_curve_empty (counter)
% % %                     % Loading data
% % %                     input_filename = sprintf('%s_%s_nc_%i.mat', gene_name, dataset_name, nuc_cyc);
% % %                     input_full_path = strcat(input_folder, input_filename);
% % %                     load (input_full_path);
% % % 
% % %                    % Preparing vertices to make the error fill
% % %                     t_vertices = new_time_mesh(1:bootstrap_only_each_frame:end);
% % %                     t_vertices = [t_vertices, flip(t_vertices)];
% % % 
% % %                     fluo_vertices = [...
% % %                         bin_fluo_data_shifted_mean(1:bootstrap_only_each_frame:end, bin_index_for_plot) + bin_fluo_data_shifted_confidence_intervals(1:bootstrap_only_each_frame:end, bin_index_for_plot,1);...
% % %                         flip(bin_fluo_data_shifted_mean(1:bootstrap_only_each_frame:end, bin_index_for_plot) - bin_fluo_data_shifted_confidence_intervals(1:bootstrap_only_each_frame:end, bin_index_for_plot, 2))]/fluo_per_polymerase;
% % % 
% % %                     %% Some of the fluo vertices can be NaN. So I need to remove those indices from all three vectors
% % %                     not_nan_indices = ~isnan(fluo_vertices);
% % %                     t_vertices = t_vertices(not_nan_indices);
% % %                     fluo_vertices = fluo_vertices(not_nan_indices);
% % % 
% % % 
% % %                     %% Creating the fill
% % % %                     fill(t_vertices, fluo_vertices', saved_curves_colors(counter, :), 'EdgeColor', 'None',...
% % % %                         'FaceAlpha', fill_opacity);
% % %                 end;
% % %                 
% % %                 counter = counter + 1;
% % %             end;
% % % 
% % % 
% % % 
% % %     %% Adjusting plot visualization parameters
% % % %     xlim(x_lim_vector);
% % %     ylim(y_lim_vector);
% % %     xlabel('Time (min)');
% % % %     ylabel(sprintf('$\\overline{N}$ (AP = %.2f)', AP_values_for_plot_array(AP_ind)), 'interpreter', 'latex');
% % %     

% % %     
% % %         end;
% % %     end;
% % % end;
% end;




% % %% Adding theoretical slope nc 14
% % figure(fig_hand_nc14);
% % % Maximal slope
% % forced_start_nc_time = 50;
% % theor_time_mesh = forced_start_nc_time:new_time_step:x_lim_vector(2);
% % theor_slope_func = @(t) k/(1+sqrt(l))^2 * t; % k in min^{-1}, t in min
% % theor_slope_data = theor_slope_func(theor_time_mesh-forced_start_nc_time);
% % plot(theor_time_mesh, theor_slope_data, '--', 'LineWidth', 1, 'color', 'black');
% %  % Slow slope
% % theor_slope_func = @(t) abortive_injection_slope_modifier * k/(1+sqrt(l))^2 * t; % k in min^{-1}, t in min
% % theor_slope_data = theor_slope_func(theor_time_mesh-forced_start_nc_time);
% % plot(theor_time_mesh, theor_slope_data, '-', 'LineWidth', 1, 'color', 'black');
% % % Max number
% % x_lim = xlim();
% % theor_max_number = L/(l+sqrt(l));
% % plot(x_lim, theor_max_number * ones(1,2), '--', 'LineWidth', 1, 'color', 'black');


% % % title ('nc 13');
% % % % NC 14
% % % figure(fig_hand_nc14);
% % % x_lim_vec = 49.5 + [0, 15];
% % % % xlim(x_lim_vec);
% % % title ('nc 14');








