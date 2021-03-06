%% Output the mean normalized values of the slope and steady-state polymerase numbers to be used in the LaTeX article.



%% Constants 
constants;
output_filename_s_NSS = 'table_s_NSS_data.tex';
output_filename_tau = 'table_tau_data.tex';
str_TASEP_norm_s_prediction = '0.34';
str_TASEP_norm_N_prediction = '0.30';



% %% Initialize
% tau



%% Open the output files
output_full_path = strcat(output_figures_folder, output_filename_s_NSS);
output_file_s_NSS = fopen(output_full_path, 'w', 'n', 'UTF-8');
output_full_path = strcat(output_figures_folder, output_filename_tau);
output_file_tau = fopen(output_full_path, 'w', 'n', 'UTF-8');



%% Calculate mean and error for the theoretical prediction of normalized Nss/NMax
% k is bp/min. Convert tau to mins
tau = tau_exp / 60;
dtau = dtau_exp / 60;

% Calculate the mean prediction
Ntheor_norm = sqrt(l) * (1 + sqrt(l)) / (k * tau + l - 1);

% Calculate the error. Remember both tau and k must be in mins
der_Ntheor_norm_l = (1 + 1 / 2 / sqrt(l) - Ntheor_norm) / (k * tau + l - 1);
der_Ntheor_norm_k = - Ntheor_norm * tau / (k * tau + l - 1);
der_Ntheor_norm_tau = - Ntheor_norm * k / (k * tau + l - 1);
d_Ntheor_norm = ((der_Ntheor_norm_l * dl)^2 + (der_Ntheor_norm_k * dk)^2 + (der_Ntheor_norm_tau * dtau)^2) ^ (1/2);

% Format strings
[str_Ntheor_norm, str_d_Ntheor_norm] = format_error_strings(Ntheor_norm, d_Ntheor_norm);



%% Calculate mean and error for the theoretical prediction of normalized s/sMax

% Calculate the mean prediction
stheor_norm = (k * tau - 1) * (1 + sqrt(l)) ^ 2 / k / tau / (k * tau + l - 1);

% Calculate the error. Remember both tau and k must be in mins
der_stheor_norm_l = ((1 + sqrt(l)) * (-1 + k * tau) * (-1 - sqrt(l) + k * tau)) / (k * sqrt(l) * tau * (-1 + l + k * tau) ^ 2);
der_stheor_norm_k = ((1 + sqrt(l)) ^ 2 * (l - (-1 + k * tau) ^ 2)) / (k ^ 2 * tau * (-1 + l + k * tau) ^ 2);
der_stheor_norm_tau = ((1 + sqrt(l)) ^ 2 * (l - (-1 + k * tau) ^ 2)) / (k * tau ^ 2 * (-1 + l + k * tau) ^ 2);

d_stheor_norm = ((der_stheor_norm_l * dl)^2 + (der_stheor_norm_k * dk)^2 + (der_stheor_norm_tau * dtau)^2) ^ (1/2);

% Format strings
[str_stheor_norm, str_d_stheor_norm] = format_error_strings(stheor_norm, d_stheor_norm);



%% Prepare the output for the steady-state polymerase number data

for nuc_cyc = nuc_cyc_array
    if nuc_cyc == 13
        str_data_line = '\\multirow{2}{\\mutlirowWidth}{$\\NSS/\\NMax$}';
        str_tau_line = '\\multirow{2}{\\mutlirowWidth}{$\\NSS/\\NMax$}';
    else
        str_data_line = '';
        str_tau_line = '';
    end;
    
    str_data_line = strcat(str_data_line, '\t&\t', num2str(nuc_cyc), '\t&\t$', str_Ntheor_norm, '\\pm', str_d_Ntheor_norm, '$');
    str_tau_line = strcat(str_tau_line, '\t&\t', num2str(nuc_cyc));
    % Cycle through genes
    for gene_ind = 1:length(gene_names_array)
        gene_name = gene_names_array{gene_ind};
        % Select gene length
        L = L_array(gene_ind);
        dL = dL_array(gene_ind);
        % Cacluate the theoretical maximum
        Ntheor = L / sqrt(l) / (1+sqrt(l));
        der_Ntheor_L = Ntheor / L;
        der_Ntheor_l = - Ntheor * (1+2*sqrt(l)) / 2 / l / (1+sqrt(l));
        d_Ntheor = ((der_Ntheor_L * dL)^2 + (der_Ntheor_l * dl)^2)^(1/2);
            
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
                    [str_N_norm, str_d_Nnorm] = format_error_strings(N_norm, d_Nnorm);
                    str_data_current = sprintf('$%s\\\\pm%s$', str_N_norm, str_d_Nnorm);
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
                    [str_tau, str_d_tau] = format_error_strings(tau, d_tau);
                    str_tau_current = sprintf('$%s\\\\pm%s$', str_tau, str_d_tau);
                    
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
% Cacluate the theoretical maximum and its error
stheor = k/(1+sqrt(l))^2;
stheor_der_k = stheor/k;
stheor_der_l = stheor / sqrt(l) / (1 + sqrt(l));
d_stheor = ((stheor_der_k * dk)^2 + (stheor_der_l * dl)^2)^(1/2);
for nuc_cyc = nuc_cyc_array
    if nuc_cyc == 13
        str_data_line = '\\multirow{2}{\\mutlirowWidth}{$s/\\sMax$}';
        str_tau_line = '\\multirow{2}{\\mutlirowWidth}{$s/\\sMax$}';
    else
        str_data_line = '';
        str_tau_line = '';
    end;
    
    str_data_line = strcat(str_data_line, '\t&\t', num2str(nuc_cyc), '\t&\t$', str_stheor_norm, '\\pm', str_d_stheor_norm, '$');
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
                % Identify non nan slopes
                not_nan_slopes = slopes_array(~isnan(slopes_array));
                % Calculate and format mean and std
                if ~isempty(not_nan_slopes)
                    s = mean(not_nan_slopes);
                    ds = std(not_nan_slopes);
                    s_norm = s / stheor;
                    d_snorm = ((ds / stheor)^2 + (s / stheor * d_stheor / stheor)^2) ^ (1/2);
                    [str_s_norm, str_d_snorm] = format_error_strings(s_norm, d_snorm);
                    str_data_current = sprintf('$%s\\\\pm%s$', str_s_norm, str_d_snorm);
                    
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
                    [str_tau, str_d_tau] = format_error_strings(tau, d_tau);
                    str_tau_current = sprintf('$%s\\\\pm%s$', str_tau, str_d_tau);
                    
                    
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



%% Close output files
fclose(output_file_s_NSS);
fclose(output_file_tau);

























