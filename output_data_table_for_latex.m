


%% This function outputs the mean normalized values of the slope and steady-state polymerase numbers to be used 
%% in the LaTeX article


%% Constants 
bin_width_slope_diml = 1;   % bin width slope
bin_width_max_number = 5;   % bin width NSS
half_width_max_rgn_ind = 2; % Number of points to consider around maximum when calculating the value of
                            % the plateau of concentration
input_folder = './processed_data/';
output_figures_folder = './figures_for_article/';
output_filename = 'table_s_NSS_data.tex';
L = 5400;      % Gene length in base pairs
k = 26*60;     % Bulk jump rate in bp/min
l = 45;        % Polymerase footprint in bp
str_TASEP_norm_s_prediction = '0.34';
str_TASEP_norm_N_prediction = '0.30';
dataset_name_array = {'bac', 'no_primary', 'no_shadow'};
gene_names_array = {'HunchBack', 'SNAIL', 'Knirps'};
nuc_cyc_array = [13, 14];


%% Old Constants                                                        
% dataset_name = 'bac';
% gene_name = 'HunchBack';
% xlim_nc_13 = [34, 54];
% xlim_nc_14 = [48, 110];
% ylim_nc_13 = [0, 107];
% ylim_nc_14 = [0, 107];
% theory_folder = '/media/aserov/DATA/Google Drive/MATLAB/13_Open_1D_TASEP_with_promoter_and_pause/histogram_data/';
% x_lim = [0, 30];
% y_lim_factor = 1.03;
% hist_color = [117, 145, 41] / 255;
% bl_save_figures = true;

% Open the output file
output_full_path = strcat(output_figures_folder, output_filename);
output_file = fopen(output_full_path, 'w', 'n', 'UTF-8');


%% Prepare the output for the steady-state polymerase number data
for nuc_cyc = nuc_cyc_array
    if nuc_cyc == 13
        str_line = '\\multirow{2}{\\mutlirowWidth}{$\\NSS/\\NMax$}';
    else
        str_line = '';
    end;
    
    str_line = strcat(str_line, '\t&\t', num2str(nuc_cyc), '\t&\t', str_TASEP_norm_N_prediction);
    % Cycle through genes
    for gene_ind = 1:length(gene_names_array)
        gene_name = gene_names_array{gene_ind};
        % Cycle through constructs
        for dataset_ind = 1:length(dataset_name_array)
            str_current = '--';
            dataset_name = dataset_name_array{dataset_ind};
            % Load data
            input_filename = sprintf('%s_%s_nc_%i.mat', gene_name, dataset_name, nuc_cyc);
            input_full_path = strcat(input_folder, input_filename);
            if exist(input_full_path, 'file')
                load (input_full_path);

                %% Calculate N_SS
                steady_state_number = calculate_N_steady_state(ms2_combined, slopes_array, mins_per_frame,...
                    forced_start_nc_time, intersct_array, init_slope_length, half_width_max_rgn_ind, fluo_per_polymerase);
                % Cacluate the theoretical maximum
                theoretical_SS_number = L / sqrt(l) / (1+sqrt(l));
                % Identify non nan slopes
                not_nan_norm_N_SS_number = steady_state_number(~isnan(steady_state_number))/theoretical_SS_number;
                % Calculate and format mean and std
                if ~isempty(not_nan_norm_N_SS_number)
                    str_current = sprintf('$%.2f\\\\pm%.2f$', mean(not_nan_norm_N_SS_number), std(not_nan_norm_N_SS_number));
                end;
            end;
            % Print out the number
            % display(str_current);
            % Add to the output string
            str_line = strcat(str_line, '\t&\t', str_current);
        end;
    end;
%     display(str_line);
    % Output the string to file
    if nuc_cyc == 13
        str_line = strcat(str_line, '\n\\\\\n');
    else
        str_line = strcat(str_line, '\n\\vspace {2mm}\n\\\\\n');
    end;
    fprintf(output_file, str_line);
end;


%% Prepare the output for the slopes data
for nuc_cyc = nuc_cyc_array
    if nuc_cyc == 13
        str_line = '\\multirow{2}{\\mutlirowWidth}{$s/\\sMax$}';
    else
        str_line = '';
    end;
    
    str_line = strcat(str_line, '\t&\t', num2str(nuc_cyc), '\t&\t', str_TASEP_norm_s_prediction);
    % Cycle through genes
    for gene_ind = 1:length(gene_names_array)
        gene_name = gene_names_array{gene_ind};
        % Cycle through constructs
        for dataset_ind = 1:length(dataset_name_array)
            str_current = '--';
            dataset_name = dataset_name_array{dataset_ind};
            % Load data
            input_filename = sprintf('%s_%s_nc_%i.mat', gene_name, dataset_name, nuc_cyc);
            input_full_path = strcat(input_folder, input_filename);
            if exist(input_full_path, 'file')
                load (input_full_path);

                %% Calculate slopes
                % Cacluate the theoretical maximum
                theoretical_slope = k/(1+sqrt(l))^2;
                % Identify non nan slopes
                not_nan_norm_slopes = normalized_slopes_array(~isnan(normalized_slopes_array));
                % Calculate and format mean and std
                if length(not_nan_norm_slopes) > 0
                    str_current = sprintf('$%.2f\\\\pm%.2f$', mean(not_nan_norm_slopes), std(not_nan_norm_slopes));
                end;
            end;
            % Print out the number
            % display(str_current);
            % Add to the output string
            str_line = strcat(str_line, '\t&\t', str_current);
        end;
    end;
    display(str_line);
    % Output the string to file
    str_line = strcat(str_line, '\n\\\\\n');
    fprintf(output_file, str_line);
end;





% Close output file
fclose(output_file);

























