


%% This function outputs the mean normalized values of the slope and steady-state polymerase numbers to be used 
%% in the LaTeX article


%% Constants 
bin_width_slope_diml = 1;   % bin width slope
bin_width_max_number = 5;   % bin width NSS
half_width_max_rgn_ind = 2; % Number of points to consider around maximum when calculating the value of
                            % the plateau of concentration
input_folder = './processed_data/';
output_figures_folder = './figures_for_article/';
output_filename_s_NSS = 'table_s_NSS_data.tex';
output_filename_tau = 'table_tau_data.tex';
L = 5400;       % Gene length in base pairs
dL = 0;         % Error in gene length
k = 26*60;      % Bulk jump rate in bp/min
dk = 2 * 60;    % Error in k
l = 45;         % Polymerase footprint in bp
dl = 5;         % Error in l
str_TASEP_norm_s_prediction = '0.34';
str_TASEP_norm_N_prediction = '0.30';
dataset_name_array = {'bac', 'no_primary', 'no_shadow'};
gene_names_array = {'HunchBack', 'SNAIL', 'Knirps'};
nuc_cyc_array = [13, 14];

% %% Private derivatives
% tau_N_der_k = 


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

% Open the output files
output_full_path = strcat(output_figures_folder, output_filename_s_NSS);
output_file_s_NSS = fopen(output_full_path, 'w', 'n', 'UTF-8');
output_full_path = strcat(output_figures_folder, output_filename_tau);
output_file_tau = fopen(output_full_path, 'w', 'n', 'UTF-8');


%% Prepare the output for the steady-state polymerase number data
 % Cacluate the theoretical maximum
Ntheor = L / sqrt(l) / (1+sqrt(l));
der_Ntheor_L = Ntheor / L;
der_Ntheor_l = - Ntheor * (1+2*sqrt(l)) / 2 / l / (1+sqrt(l));
d_Ntheor = ((der_Ntheor_L * dL)^2 + (der_Ntheor_l * dl)^2)^(1/2);

for nuc_cyc = nuc_cyc_array
    if nuc_cyc == 13
        str_data_line = '\\multirow{2}{\\mutlirowWidth}{$\\NSS/\\NMax$}';
        str_tau_line = '\\multirow{2}{\\mutlirowWidth}{$\\NSS/\\NMax$}';
    else
        str_data_line = '';
        str_tau_line = '';
    end;
    
    str_data_line = strcat(str_data_line, '\t&\t', num2str(nuc_cyc), '\t&\t', str_TASEP_norm_N_prediction);
    str_tau_line = strcat(str_tau_line, '\t&\t', num2str(nuc_cyc));
    % Cycle through genes
    for gene_ind = 1:length(gene_names_array)
        gene_name = gene_names_array{gene_ind};
        % Cycle through constructs
        for dataset_ind = 1:length(dataset_name_array)
            str_data_current = '--';
            str_tau_current = '--';
            dataset_name = dataset_name_array{dataset_ind};
            % Load data
            input_filename = sprintf('%s_%s_nc_%i.mat', gene_name, dataset_name, nuc_cyc);
            input_full_path = strcat(input_folder, input_filename);
            if exist(input_full_path, 'file')
                load (input_full_path);
                %% Calculate N_SS
                steady_state_number = calculate_N_steady_state(ms2_combined, slopes_array, mins_per_frame,...
                    forced_start_nc_time, intersct_array, init_slope_length, half_width_max_rgn_ind, fluo_per_polymerase);
                % Identify non nan slopes
                not_nan_N_SS_number = steady_state_number(~isnan(steady_state_number));
                % Calculate and format mean and std
                if ~isempty(not_nan_N_SS_number)
                    N = mean(not_nan_N_SS_number);
                    dN = std(not_nan_N_SS_number);
                    %
                    N_norm = N/Ntheor;
                    d_Nnorm = ((dN/Ntheor)^2 + (N/Ntheor * d_Ntheor / Ntheor)^2) ^ (1/2);
                    str_data_current = sprintf('$%.2f\\\\pm%.2f$', N_norm, d_Nnorm);
                    % Calculate tau
                    tau = (sqrt(l) * (1 + sqrt(l)) * Ntheor / N - l + 1)/k;         
                    % Calculate the error in tau estimate
                    tau_N_der_k = - tau/k;
                    tau_N_der_Ntheor = sqrt(l) * (1 + sqrt(l)) / N / k;
                    tau_N_der_N = - sqrt(l) * (1 + sqrt(l)) * Ntheor / N^2/k;
                    tau_N_der_l = (1/2/sqrt(l) * Ntheor / N + Ntheor/N - 1)/k;
                    %
                    d_tau = ((tau_N_der_k * dk)^2 + (tau_N_der_Ntheor * d_Ntheor)^2 + (tau_N_der_N * dN)^2 + (tau_N_der_l * dl)^2) ^ (1/2);
                    % Convert to seconds
                    tau = tau * 60;
                    d_tau = d_tau * 60;
                    str_tau_current = sprintf('$%.1f\\\\pm%.1f$', tau, d_tau);
                    
                end;
            end;
            % Print out the number
            % display(str_current);
            % Add to the output string
            str_data_line = strcat(str_data_line, '\t&\t', str_data_current);
            str_tau_line = strcat(str_tau_line, '\t&\t', str_tau_current);
        end;
    end;
%     display(str_line);
    % Output the string to file
    if nuc_cyc == 13
        str_data_line = strcat(str_data_line, '\n\\\\\n');
        str_tau_line = strcat(str_tau_line, '\n\\\\\n');
    else
        str_data_line = strcat(str_data_line, '\n\\vspace {2mm}\n\\\\\n');
        str_tau_line = strcat(str_tau_line, '\n\\vspace {2mm}\n\\\\\n');
    end;
    fprintf(output_file_s_NSS, str_data_line);
    fprintf(output_file_tau, str_tau_line);
end;


%% Prepare the output for the slopes data
for nuc_cyc = nuc_cyc_array
    if nuc_cyc == 13
        str_data_line = '\\multirow{2}{\\mutlirowWidth}{$s/\\sMax$}';
        str_tau_line = '\\multirow{2}{\\mutlirowWidth}{$s/\\sMax$}';
    else
        str_data_line = '';
        str_tau_line = '';
    end;
    
    str_data_line = strcat(str_data_line, '\t&\t', num2str(nuc_cyc), '\t&\t', str_TASEP_norm_s_prediction);
    str_tau_line = strcat(str_tau_line, '\t&\t', num2str(nuc_cyc));
    % Cycle through genes
    for gene_ind = 1:length(gene_names_array)
        gene_name = gene_names_array{gene_ind};
        % Cycle through constructs
        for dataset_ind = 1:length(dataset_name_array)
            str_data_current = '--';
            dataset_name = dataset_name_array{dataset_ind};
            % Load data
            input_filename = sprintf('%s_%s_nc_%i.mat', gene_name, dataset_name, nuc_cyc);
            input_full_path = strcat(input_folder, input_filename);
            if exist(input_full_path, 'file')
                load (input_full_path);

                %% Calculate slopes
                % Cacluate the theoretical maximum and its error
                stheor = k/(1+sqrt(l))^2;
                stheor_der_k = s/k;
                stheor_der_l = stheor / sqrt(l) / (1 + sqrt(l));
                d_stheor = ((stheor_der_k * dk)^2 + (stheor_der_l * dl)^2)^(1/2);
                % Identify non nan slopes
                not_nan_slopes = slopes_array(~isnan(slopes_array));
                % Calculate and format mean and std
                if ~isempty(not_nan_slopes)
                    s = mean(not_nan_slopes);
                    ds = std(not_nan_slopes);
                    str_data_current = sprintf('$%.2f\\\\pm%.2f$', s/stheor, ds/stheor);
                    
                    % Calculate tau
                    tau = (1/(2*k^2*s)) * (k*(s - l*s + (1 + sqrt(l))^2 * stheor) + stheor *...
                        sqrt((k^2*(1 + sqrt(l))^2 * ((-1 + sqrt(l))^2*s^2 - 2*(1 + l)*s * stheor + (1 + sqrt(l))^2 *...
                        stheor^2))/stheor^2));
                    % Calculate the error in tau estimate
                    % Need derivatives
                    tau_s_der_k = - tau/k;
                    tau_s_der_l = -((s - stheor - stheor/sqrt(l))/(2*k*s)) + ((1 + sqrt(l))*(-s + stheor) * ...
                                (stheor + l*(-s + stheor) + sqrt(l)*(s + 2*stheor))) / (2*sqrt(l)*s*stheor*sqrt((1/stheor^2) * ...
                                (k^2*(1 + sqrt(l))^2*((-1 + sqrt(l))^2 * s^2 - 2*(1 + l)*s*stheor + (1 + sqrt(l))^2*stheor^2))));
                    tau_s_der_s = -(((1 + sqrt(l))^2*((-k)*(1 + l)*s +k*(1 + sqrt(l))^2*stheor + stheor*sqrt((1/stheor^2)*...
                                (k^2*(1 + sqrt(l))^2*((-1 + sqrt(l))^2*s^2 - 2*(1 + l)*s*stheor + (1 + sqrt(l))^2*...
                                stheor^2)))))/(2*k*s^2*sqrt((1/stheor^2)*(k^2*(1 + sqrt(l))^2*((-1 + sqrt(l))^2*s^2 - 2*(1 + l)*s*...
                                stheor + (1 + sqrt(l))^2*stheor^2)))));
                    tau_s_der_stheor = ((1 + sqrt(l))^2*((-k)*(1 + l)*s + k*(1 + sqrt(l))^2*stheor + stheor*sqrt((1/stheor^2)*...
                                (k^2*(1 + sqrt(l))^2*((-1 + sqrt(l))^2*s^2 - 2*(1 + l)*s*stheor + (1 + sqrt(l))^2*stheor^2)))))/...
                                (2*k*s*stheor*sqrt((1/stheor^2)*(k^2*(1 + sqrt(l))^2*((-1 + sqrt(l))^2*s^2 - 2*(1 + l)*s*...
                                stheor + (1 + sqrt(l))^2*stheor^2))));
                    %
                    d_tau = ((tau_s_der_k * dk)^2 + (tau_s_der_l * dl)^2 + (tau_s_der_s * ds)^2 + (tau_s_der_stheor * d_stheor)^2) ^ (1/2);
                    
                    % Convert to seconds
                    tau = tau * 60;
                    d_tau = d_tau * 60;
                    str_tau_current = sprintf('$%.1f\\\\pm%.1f$', tau, d_tau);
                    
                    
                end;
            end;
            % Print out the number
            % display(str_current);
            % Add to the output string
            str_data_line = strcat(str_data_line, '\t&\t', str_data_current);
            str_tau_line = strcat(str_tau_line, '\t&\t', str_tau_current);
        end;
    end;
%     display(str_data_line);
    % Output the string to file
    str_data_line = strcat(str_data_line, '\n\\\\\n');
    fprintf(output_file_s_NSS, str_data_line);
    %
    str_tau_line = strcat(str_tau_line, '\n\\\\\n');
    fprintf(output_file_tau, str_tau_line);
end;





% Close output files
fclose(output_file_s_NSS);
fclose(output_file_tau);

























