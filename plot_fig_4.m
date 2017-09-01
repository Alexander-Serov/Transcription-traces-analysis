



%% Plotting the maximal fluorescence curve for each gene





%% Constants
dataset_names_array = {'bac', 'no_primary', 'no_shadow'};
input_folder = './processed_data/';
fill_opacity = 0.2;
x_lim_vector = [48, 80];
AP_values_for_plot_array = 0.15:0.05:0.55;
% bl_save_figures = true;


% Initializing the figure
fig_hand = figure(14);
% set_my_fig_size (fig_hand);
clf;
% hold on;


for AP_ind = 1:length(AP_values_for_plot_array)
    subplot(3, ceil(length(AP_values_for_plot_array)/3), AP_ind);
    hold on;

    AP_value_for_plot = AP_values_for_plot_array(AP_ind);



    %% Initialization
    bin_index_for_plot = ceil((AP_value_for_plot-bins_borders(1))/bin_width);



    %% Parsing different genes first to plot the mean values
    for construct_index = 1:length(dataset_names_array)
        dataset_name = dataset_names_array{construct_index};


        % Loading data
        input_filename = sprintf('%s_%s_nc_%i.mat', gene_name, dataset_name, nuc_cyc);
        input_full_path = strcat(input_folder, input_filename);
        load (input_full_path);

    % % %     % Determining the maximal fluorescence bin for each construct
    % % %     [~, max_bin_index] =  max(max(bin_fluo_data_shifted_mean, [], 1));

        % Plotting the bin we found
        plot(new_time_mesh, bin_fluo_data_shifted_mean(:, bin_index_for_plot)/fluo_per_polymerase,...
                'color', my_constructs_color_sequence(construct_index, :), 'LineWidth', 2);

    end;


    % Adding the legend
    legend_string = cell(length(dataset_names_array), 1);
    for i = 1:length(dataset_names_array)
        legend_string{i} = sprintf('%s, AP = %.2f', dataset_names_array{i}, bins_centers(bin_index_for_plot));
        legend_string{i} = strrep(legend_string{i}, '_', ' ');
    end;
%     if AP_ind == 1
        legend(legend_string, 'Location', 'northeast', 'FontSize', 12);
%     end;


    %% Now parsing again to plot the fill.
    % Two for loops are required for correct legends

    for construct_index = 1:length(dataset_names_array)
        dataset_name = dataset_names_array{construct_index};


        % Loading data
        input_filename = sprintf('%s_%s_nc_%i.mat', gene_name, dataset_name, nuc_cyc);
        input_full_path = strcat(input_folder, input_filename);
        load (input_full_path);

    % % %     % Determining the maximal fluorescence bin for each construct
    % % %     [~, max_bin_index] =  max(max(bin_fluo_data_shifted_mean, [], 1));

        % Preparing vertices to make the error fill
        t_vertices = new_time_mesh(1:bootstrap_only_each_frame:end);
        t_vertices = [t_vertices, flip(t_vertices)];

        fluo_vertices = [...
            bin_fluo_data_shifted_mean(1:bootstrap_only_each_frame:end, bin_index_for_plot) + bin_fluo_data_shifted_confidence_intervals(1:bootstrap_only_each_frame:end, bin_index_for_plot,1);...
            flip(bin_fluo_data_shifted_mean(1:bootstrap_only_each_frame:end, bin_index_for_plot) - bin_fluo_data_shifted_confidence_intervals(1:bootstrap_only_each_frame:end, bin_index_for_plot, 2))]/fluo_per_polymerase;

        %% Some of the fluo vertices can be NaN. So I need to remove those indices from all three vectors
        not_nan_indices = ~isnan(fluo_vertices);
        t_vertices = t_vertices(not_nan_indices);
        fluo_vertices = fluo_vertices(not_nan_indices);


        %% Creating the fill
        fill(t_vertices, fluo_vertices', my_constructs_color_sequence(construct_index, :), 'EdgeColor', 'None',...
            'FaceAlpha', fill_opacity);

    end;


    %% Adding maximal injection prediciton

    % The 'expected l' curve
    theor_time_mesh = forced_start_nc_time:new_time_step:x_lim_vector(2);
    theor_slope_func = @(t) k/(1+sqrt(l))^2 * t; % k in min^{-1}, t in min
    theor_slope_data = theor_slope_func(theor_time_mesh-forced_start_nc_time);
    plot(theor_time_mesh, theor_slope_data, '--', 'LineWidth', 2, 'color', 'black');

    % Adding the fill
    x_vertices = forced_start_nc_time + [0, 20, 20, 0];
    y_vertices = [0,...
        k/(1+sqrt(l_min))^2 * (x_vertices(2) - x_vertices(1)),...
        k/(1+sqrt(l_max))^2 * (x_vertices(3) - x_vertices(1)),...
        0];
    fill(x_vertices, y_vertices, 'k', 'FaceAlpha', fill_opacity, 'EdgeColor', 'None');


    %% Adding maximal particle number prediction
    plot(x_lim_vector, ones(1,2) * L/(l+sqrt(l)), '-', 'color', 'black', 'LineWidth', 2);



    %% Adjusting plot visualization parameters

    xlim(xlim_vector);
    ylim([0, 80]);

    xlabel('Time (min)');
    ylabel('$\overline{N_p}$', 'interpreter', 'latex');


    hold off;


end;



% Saving the figure
output_filename_svg = sprintf('04_%s_Max_luminescence_dif_constructs_nc_%i.png', gene_name, nuc_cyc);
output_filename_pdf = sprintf('04_%s_Max_luminescence_dif_constructs_nc_%i.pdf', gene_name, nuc_cyc);
output_full_path_svg = strcat(output_figures_folder, output_filename_svg);
output_full_path_pdf = strcat(output_figures_folder, output_filename_pdf);
if bl_save_figures
    fig_hand.PaperPositionMode = 'auto';
    saveas(fig_hand, output_full_path_svg, 'png');
%     saveas(fig_hand, output_full_path_pdf, 'pdf');
    display('Figure saved');
end;



