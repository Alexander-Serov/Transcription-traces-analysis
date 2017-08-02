



%% Constants
dataset_names_array = {'bac', 'no_primary', 'no_shadow'};
dataset_short_names_array = {'bac', 'no_pr', 'no_sh'};
gene_names_array = {'HunchBack', 'Knirps', 'SNAIL'};
gene_short_names_array = {'Hb', 'Kn', 'Sn'};
input_folder = './processed_data/';
fill_opacity = 0.13;
% abortive_injection_slope_modifier = 0.34;
% y_lim_vector = 'auto'; % [0, 1];
y_lim_vector = [0, 1];
% y_lim_vector = [0.08, 0.62];
bl_legend_on = true;
% bl_save_figures = false;
% bl_save_figures = true;
% image_size = [1000, 300];
% fig_size = [360, 370];
fig_side = 400;
AP_rand_shift_size = 0.006;



%% Initialization
max_curve_number = 2 * length(dataset_names_array);
bl_curve_empty = zeros(length(gene_names_array), max_curve_number);
saved_curves_colors = zeros(length(gene_names_array), max_curve_number, 3);
rng('shuffle');

%% Initializing the figure
fig_hand = figure(8);
set_my_fig_size (fig_hand, fig_side);
clf;
% hold on;


%% Parsing different constructs
shifts = 2* AP_rand_shift_size*(rand(1, 2*length(dataset_names_array)) - 1/2);

for gene_index = 1:1% length(gene_names_array)
%     subplot(1, 3, gene_index);
    hold on;
    % Proceeding to plot
    legend_counter = 1;
    counter = 1;
    legend_string = cell(max_curve_number, 1);
    gene_name = gene_names_array{gene_index};
    min_AP = NaN;
    max_AP = NaN;
    for nuc_cyc = [13, 14]
        for construct_index = 1:length(dataset_names_array)
            dataset_name = dataset_names_array{construct_index};
            % Loading data
            input_filename = sprintf('%s_%s_nc_%i.mat', gene_name, dataset_name, nuc_cyc);
            input_full_path = strcat(input_folder, input_filename);
            load (input_full_path);
            number_of_not_nan_points = sum(~isnan(bin_normalized_slopes_mean));
            % Loading error bars
            error_bar = bin_normalized_slopes_confidence_intervals(:, 1:2);
            error_bar(:, 1) = error_bar(:, 1) - bin_normalized_slopes_mean';
            error_bar(:, 2) = error_bar(:, 2) - bin_normalized_slopes_mean';
            if number_of_not_nan_points > 0
                h_curve = errorbar(bins_centers + shifts(legend_counter),...
                    bin_normalized_slopes_mean,...
                        error_bar(:, 1), error_bar(:, 2), '-o',...
                        ... %'color', my_constructs_color_sequence(construct_index, :), 
                        'LineWidth', 2, 'MarkerSize', 8);
%                 h_curve = plot(bins_centers, bin_normalized_slopes_mean,...
%                         ... %'color', my_constructs_color_sequence(construct_index, :), 
%                         'LineWidth', 3);
                saved_curves_colors(gene_index, counter, :) = get(h_curve, 'Color');
                % Detecting the plot x limits
                valid_AP_values = bins_centers(~isnan(bin_normalized_slopes_mean));
                min_AP = min(min_AP, min(valid_AP_values));
                max_AP = max(max_AP, max(valid_AP_values));
                

                %% Preparing the legend
                legend_string{legend_counter} = sprintf('nc %i, %s', nuc_cyc, dataset_short_names_array{construct_index});
                legend_string{legend_counter} = strrep(legend_string{legend_counter}, '_', ' ');
                legend_counter = legend_counter + 1;
            else
                bl_curve_empty(gene_index, counter) = 1;
            end;
            counter = counter + 1;

        end;
    end;
    
    %% Adding the legend
    % Removing empty legends
    legend_string(cellfun(@isempty,legend_string)) = [];
    if bl_legend_on
        legend(legend_string, 'FontSize', 12, 'location', 'best');
    end;
end;    % gene_index


%% Adjusting the plot
%% Plot the slow sites prediction
x_lim_vector = xlim() + 0.05 .* [-1, 1];
theor_slope = (k/60*tau_exp - 1) * (1 + sqrt(l))^2 / (k/60 * tau_exp) / (k/60 * tau_exp + l - 1);     % this slope is already normalized to k and converted to mins
plot(x_lim_vector, theor_slope * ones(2,1), '-', 'LineWidth', 2, 'color', 'black');
%% Adjusting the plotting
% xlim(x_lim_vector);
ylim(y_lim_vector);
xlabel('AP position');
ylabel('s/s_{max}');
title(gene_short_names_array{gene_index});


% % % % % %% Now parsing again to plot the fill.
% % % % % % Two for loops are required for correct legends
% % for gene_index = 1:1; length(gene_names_array)
% % %     subplot(1, 3, gene_index);
% %     hold on;
% %     % Proceeding to plot
% %     legend_counter = 1;
% %     counter = 1;
% %     legend_string = cell(max_curve_number, 1);
% %     gene_name = gene_names_array{gene_index};
% %     min_AP = NaN;
% %     max_AP = NaN;
% %     for nuc_cyc = [13, 14]
% %         for construct_index = 1:length(dataset_names_array)
% %             dataset_name = dataset_names_array{construct_index};
% %             if ~bl_curve_empty(gene_index, counter)
% %                 % Loading data
% %                 input_filename = sprintf('%s_%s_nc_%i.mat', gene_name, dataset_name, nuc_cyc);
% %                 input_full_path = strcat(input_folder, input_filename);
% %                 load (input_full_path);
% %                 % Preparing vertices to make the error fill
% %                 AP_vertices = bins_centers;
% %                 AP_vertices = [AP_vertices, flip(AP_vertices)];
% % 
% %                 slope_vertices = [...
% %                     bin_normalized_slopes_confidence_intervals(:, 1);...
% %                     flip(bin_normalized_slopes_confidence_intervals(:, 2))];
% % 
% %                 %% Some of the fluo vertices can be NaN. So I need to remove those indices from all three vectors
% %                 not_nan_indices = ~isnan(slope_vertices) & slope_vertices > 0;
% %                 AP_vertices = AP_vertices(not_nan_indices);
% %                 slope_vertices = slope_vertices(not_nan_indices);
% % 
% % 
% %                 %% Creating the fill
% %                 fill(AP_vertices, slope_vertices', saved_curves_colors(gene_index, counter, :), 'EdgeColor', 'None',...
% %                     'FaceAlpha', fill_opacity);
% %             end;
% %             counter = counter + 1;
% %         end;
% %     end;
% % end;    % gene_index


%% Adjusting appearance
grid on;
grid minor;
box on;
% % xlim('auto');
xlim([0.14, 0.51]);
ylim([0,0.74]);





%% Saving the figure
if bl_save_figures
    output_filename_png = sprintf('03_Slope_vs_AP_article.png');
    output_full_path_png = strcat(output_figures_folder, output_filename_png);
    display('Saving figure');
    fig_hand.PaperPositionMode = 'auto';
    saveas(fig_hand, output_full_path_png, 'png');
    % PDF
    output_filename = sprintf('03_Slope_vs_AP_article.pdf');
    output_full_path = strcat(output_figures_folder, output_filename);
    figure(fig_hand);
    % Adjusting size
    set(fig_hand,'Units','Inches');
    pos = get(fig_hand,'Position');
    set(fig_hand,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);

    print(fig_hand, output_full_path, '-dpdf', '-r0');

    display('Figure saved');
end;








