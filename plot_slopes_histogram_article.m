

%% Plotting histograms of slope distributions for different genes and constructs
% Remember my saved slopes are in pol/min


%% Constants 
input_folder = './processed_data/';
theory_folder = '/Users/alexander_serov/Google Drive/MATLAB/13_Open_1D_TASEP_with_promoter_and_pause/histogram_data/';
% output_figures_folder = './figures_for_article/';
dataset_name = 'bac';
gene_name = 'HunchBack';
xlim_nc_13 = [34, 54];
xlim_nc_14 = [48, 110];
ylim_nc_13 = [0, 107];
ylim_nc_14 = [0, 107];
bin_width_slope_diml = 1;
bin_width_max_number = 5;
x_lim = [0, 30];
y_lim_factor = 1.03;
half_width_max_rgn_ind = 2; % Number of points to consider around maximum when calculating the value of
                            % the plateau of concentration
hist_color = [117, 145, 41] / 255;
% bl_save_figures = true;


%% Initialization
fprintf('\n');
std_data_length = 2;
experimental_std_data = struct;
experimental_std_data.('Nuc_Cyc') = zeros(1, std_data_length);
experimental_std_data.('Gene') = cell(1, std_data_length);
experimental_std_data.('Construct') = cell(1, std_data_length);
experimental_std_data.('STD_Slope') = zeros(1, std_data_length);
experimental_std_data.('STD_SS_Number') = zeros(1, std_data_length);


h_fig_slopes = cell(1,2);
h_fig_max_number = cell(1,2);
h_fig_slopes{1} = figure (8);
clf('reset');
h_fig_slopes{2} = figure (9);
clf('reset');
h_fig_max_number{1} = figure (10);
clf('reset');
h_fig_max_number{2} = figure (11);
clf('reset');

fig_side_ratio = 1.8/2;
fig_side = 320;
set_my_fig_size(h_fig_slopes{1}, fig_side * [fig_side_ratio, 1]);
set_my_fig_size(h_fig_slopes{2}, fig_side * [fig_side_ratio, 1]);
set_my_fig_size(h_fig_max_number{1}, fig_side * [fig_side_ratio, 1]);
set_my_fig_size(h_fig_max_number{2}, fig_side * [fig_side_ratio, 1]);


%% Plotting two nuclear cycles
std_data_counter = 1;
nuc_cyc_array = [13, 14];
for counter = 1:length(nuc_cyc_array)
    nuc_cyc = nuc_cyc_array(counter);
    % Loading data
    input_filename = sprintf('%s_%s_nc_%i.mat', gene_name, dataset_name, nuc_cyc);
    input_full_path = strcat(input_folder, input_filename);
    load (input_full_path);
    %% Plotting the injection rate slopes 
    figure(h_fig_slopes{counter});
% %     subplot(1, length(nuc_cyc_array), counter);
    hold on;
    % for i=1:ms2_spots_count
    %     if ms2_combined(i).nuc_cyc ~= nuc_cyc; continue; end;
    %     plot (ms2_combined(i).Frame * mins_per_frame, ms2_combined(i).Fluo/fluo_per_polymerase);
    % end;
    max_slope = 1.1 * max(slopes_array);
    bins_slope = 0:bin_width_slope_diml:max_slope;
    hist_hand = histogram(slopes_array, bins_slope, 'Normalization', 'pdf', 'FaceColor', hist_color);
    x_lim = [0, max_slope];
    xlim(x_lim);
    y_lim = [0, max(hist_hand.Values)];
%     y_lim = y_lim_factor * y_lim;
%     ylim(y_lim);
    xlabel('Initial injection rate s, pol/min');
    ylabel('PDF');
    title(sprintf('nc %i', nuc_cyc));
    % Adding the theoretical maximum
    theoretical_slope = k/(1+sqrt(l))^2;
    y_lims = ylim();
    plot(ones(1,2) * theoretical_slope, y_lims, '--k', 'LineWidth', 2);

        
    %% Printing out the mean value
    not_nan_slopes = slopes_array(~isnan(slopes_array));
    fprintf('Mean polymerase loading rate (nc %i): %.0f +- %.0f\n', nuc_cyc, mean(not_nan_slopes),...
        std(not_nan_slopes));
    fprintf('Mean normalized polymerase loading rate (nc %i): %.2f +- %.2f\n',...
        nuc_cyc, mean(not_nan_slopes/theoretical_slope),...
        std(not_nan_slopes/theoretical_slope));
    
% % %     %% Adding the calculated theoretical histogram
% % %     % Loading data
% % %     str_type = 'slope';
% % %     file_name = sprintf('theory_nc_%i_type_%s.mat', nuc_cyc, str_type);
% % %     full_path = strcat(input_theory_folder, file_name);
% % %     if exist(full_path, 'file')
% % %         load(full_path, '-mat', 'hist_bins_centers', 'hist_values');
% % %         % Plotting the theoretical curve
% % %         plot(hist_bins_centers, hist_values, 'k', 'LineWidth', 2);
% % %     end
    hold off;
    
    
    
    %% ===
    %% And now plotting the maximal polymerase numbers
    % This code is modified to give a better estimate of the plateau value.
    % Now we will take not just the peak value, but also 2 points on each
    % side of it and calculate the average
    ms2_count = size(ms2_combined, 2);
    steady_state_number = zeros(ms2_count, 1);
    for ms2_ind=1:ms2_count
        % Taking only traces for which we have a well defined slope
        if ~isnan(slopes_array(ms2_ind))
            % Finding index max. value
            cur_trace_length = length(ms2_combined(ms2_ind).Fluo);
            max_ind = find(ms2_combined(ms2_ind).Fluo == max(ms2_combined(ms2_ind).Fluo));
            % Determining plateau start index
            plateau_rgn_frame_indices = ((ms2_combined(ms2_ind).Frame * mins_per_frame) + (forced_start_nc_time - intersct_array(ms2_ind))) - ...
                (forced_start_nc_time + init_slope_length) > 0;
            plateau_start_index = find(plateau_rgn_frame_indices, 1);
            % Chossing a maximum of 5 indices around the max value to include into the mean plateau
            % calculation
            if plateau_start_index > 0
                % Using the plateau value only if the plateau region was detected
                rgn_start_index = max([max_ind - half_width_max_rgn_ind, plateau_start_index]);
                rgn_end_index = min([max_ind + half_width_max_rgn_ind, cur_trace_length]);
                % Calculating the plateau estimate
                cur_data = ms2_combined(ms2_ind).Fluo;
                cur_data = cur_data(rgn_start_index:rgn_end_index)/fluo_per_polymerase;
                steady_state_number(ms2_ind) = mean(cur_data);
            else
                % Ignoring data if the plateau region wasn't detected
                steady_state_number(ms2_ind) = NaN;
            end;
        else
            % If the slope is not defined, saving as NaN
            steady_state_number(ms2_ind) = NaN;
        end;
    end;
    max_max_number = 1.1 * max(steady_state_number);
    bins_max_numbers = 0:bin_width_max_number:max_max_number;
    % Plotting the histogram
    figure(h_fig_max_number{counter});
% %     subplot(1, length(nuc_cyc_array), counter);
    hold on;
    hist_hand = histogram(steady_state_number, bins_max_numbers, 'Normalization', 'pdf', 'FaceColor', hist_color);
    % Adjusting the plot appearence
    theoretical_max_max = L/(l+sqrt(l));
    x_lim = [0, 1.02 * max(max_max_number, theoretical_max_max)];
    xlim(x_lim);
    xlabel('Max. polymerase number N_{ss}');
    ylabel('PDF');
    title(sprintf('nc %i', nuc_cyc)); 
    y_lim = [0, max(hist_hand.Values)];
    y_lim = y_lim_factor * y_lim;
    ylim(y_lim);
    
    % Adding the theoretical maximum
    y_lims = ylim();
    plot(ones(1,2) * theoretical_max_max, y_lims, '--k', 'LineWidth', 2);
    hold off;
    
    %% Outputting interesting data
    fprintf('Total number of curves in bins (nc %i): %i (slope), %i (steady state)\n', nuc_cyc,...
        sum(~isnan(slopes_array)), sum(~isnan(steady_state_number)));

    %% Preparing to share the STD experimental data 
    experimental_std_data.('Nuc_Cyc')(std_data_counter) = nuc_cyc;
    experimental_std_data.('Gene'){std_data_counter} = gene_name;
    experimental_std_data.('Construct'){std_data_counter} = dataset_name;
    experimental_std_data.('STD_Slope')(std_data_counter) = std(not_nan_slopes);
    experimental_std_data.('STD_Slope_Norm')(std_data_counter) = std(not_nan_slopes/theoretical_slope);
    experimental_std_data.('STD_SS_Number')(std_data_counter) = std(steady_state_number(~isnan(steady_state_number)));
    
    % Increasing the counter
    std_data_counter = std_data_counter + 1;
    
    %% Printing out the mean value
    fprintf('Steady state polymerase number (nc %i): %.0f +- %.0f\n', nuc_cyc,...
        mean(steady_state_number(~isnan(steady_state_number))),...
        std(steady_state_number(~isnan(steady_state_number))));
    

    
end;

%% Saving the experimental STD data
display('Saving experimental STD data...');
file_name = sprintf('std_experimental_data_gene_%s_%s.mat', gene_name, dataset_name);
full_path = strcat(theory_folder, file_name);
save(full_path, 'experimental_std_data', '-mat');
display('Experimental STD data saved successfully!');






%% Saving the figures
if bl_save_figures
    display('Saving figures...');
    for counter = 1:2
        h_fig = h_fig_slopes{counter};
        % Slopes
        % Adjusting size
        h_fig.PaperPositionMode = 'auto';
        set(h_fig,'Units','Inches');
        pos = get(h_fig,'Position');
        set(h_fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);

        % PNG
        output_filename = sprintf('01a_%i_%s_slope_histograms_%s.png', counter, gene_name, dataset_name);
        output_full_path = strcat(output_figures_folder, output_filename);
        saveas(h_fig, output_full_path, 'png');

        % PDF
        output_filename = sprintf('01a_%i_%s_slope_histograms_%s.pdf', counter, gene_name, dataset_name);
        output_full_path = strcat(output_figures_folder, output_filename);


        print(h_fig, output_full_path, '-dpdf', '-r0');
    end;
    display('Figure saved');
    % Max numbers
    for counter = 1:2
        h_fig = h_fig_max_number{counter};
        % Adjusting size
        set(h_fig,'Units','Inches');
        pos = get(h_fig,'Position');
        set(h_fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
        % PNG
        output_filename = sprintf('01b_%i_%s_max_numbers_histograms_%s.png', counter, gene_name, dataset_name);
        output_full_path = strcat(output_figures_folder, output_filename);
        h_fig.PaperPositionMode = 'auto';
        saveas(h_fig, output_full_path, 'png');
        % PDF
        output_filename = sprintf('01b_%i_%s_max_numbers_histograms_%s.pdf', counter, gene_name, dataset_name);
        output_full_path = strcat(output_figures_folder, output_filename);
        figure(h_fig);


        print(h_fig, output_full_path, '-dpdf', '-r0');
    end;
    display('Figure saved');
end;







