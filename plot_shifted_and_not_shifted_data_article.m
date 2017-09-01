

%% Combining plots of shifted and not shifted data for n13 and nc14 for
% Hunchback only, bac construct


%% Constants 
input_folder = './processed_data/';
dataset_name = 'bac';
gene_name = 'HunchBack';
xlim_nc_13 = [34, 54];
xlim_nc_14 = [48, 110];
ylim_nc_13 = [0, 107];
ylim_nc_14 = [0, 107];



h_fig = figure (6);
clf;

set_my_fig_size(h_fig);

%% nc13 raw data
nuc_cyc=13;
% Loading data
input_filename = sprintf('%s_%s_nc_%i.mat', gene_name, dataset_name, nuc_cyc);
input_full_path = strcat(input_folder, input_filename);
load (input_full_path);
% Plotting
subplot(2, 2, 1);
hold on;
for i=1:ms2_spots_count
    if ms2_combined(i).nuc_cyc ~= nuc_cyc; continue; end;
    plot (ms2_combined(i).Frame * mins_per_frame, ms2_combined(i).Fluo/fluo_per_polymerase);
end;
%% Adding the theoretical prediction
% hold on;
theor_time_mesh = forced_start_nc_time:0.05:(forced_start_nc_time+5);
theor_slope_func = @(t) k/(1+sqrt(l))^2 * t; % k in min^{-1}, t in min
theor_slope_data = theor_slope_func(theor_time_mesh-forced_start_nc_time);
plot(theor_time_mesh, theor_slope_data, '-', 'LineWidth', 2, 'color', 'black');
%% Adding the maximal particle number
plot(xlim_nc_13, ones(1,2) * L/(l+sqrt(l)), ':', 'color', 'black',...
    'LineWidth', 2);
%% Adjusting the plot
xlim(xlim_nc_13);
ylim(ylim_nc_13);
xlabel('Time (min)');
ylabel('N_{pol}');
title(sprintf('Raw data, nc %i', nuc_cyc));


%% nc13 synced data
nuc_cyc=13;
% Loading data
input_filename = sprintf('%s_%s_nc_%i.mat', gene_name, dataset_name, nuc_cyc);
input_full_path = strcat(input_folder, input_filename);
load (input_full_path);
% Plotting
subplot(2, 2, 2);
hold on;
for i=1:ms2_spots_count
    if ms2_combined(i).nuc_cyc ~= nuc_cyc; continue; end;
    % Taking only those ms2 for which the calculated intersect was > 0
    if intersct_array(i)>0
        plot ((ms2_combined(i).Frame * mins_per_frame) + (forced_start_nc_time - intersct_array(i)), ...
            ms2_combined(i).Fluo/fluo_per_polymerase);        
    end;
end;
% Adding the theoretical prediction
theor_time_mesh = forced_start_nc_time:0.05:(forced_start_nc_time+5);
theor_slope_func = @(t) k/(1+sqrt(l))^2 * t; % k in min^{-1}, t in min
theor_slope_data = theor_slope_func(theor_time_mesh-forced_start_nc_time);
plot(theor_time_mesh, theor_slope_data, '-', 'LineWidth', 2, 'color', 'black');
%% Adding the maximal particle number
plot(xlim_nc_13, ones(1,2) * L/(l+sqrt(l)), ':', 'color', 'black',...
    'LineWidth', 2);
%% Adjusting the plot
xlim(xlim_nc_13);
ylim(ylim_nc_13);
xlabel('Time (min)');
ylabel('N_{pol}');
title(sprintf('Synced data, nc %i', nuc_cyc));



%% nc14 raw data
nuc_cyc=14;
% Loading data
input_filename = sprintf('%s_%s_nc_%i.mat', gene_name, dataset_name, nuc_cyc);
input_full_path = strcat(input_folder, input_filename);
load (input_full_path);
% Plotting
subplot(2, 2, 3);
hold on;
for i=1:ms2_spots_count
    if ms2_combined(i).nuc_cyc ~= nuc_cyc; continue; end;
    plot (ms2_combined(i).Frame * mins_per_frame, ms2_combined(i).Fluo/fluo_per_polymerase);
end;
%% Adding the theoretical prediction
% hold on;
theor_time_mesh = forced_start_nc_time:0.05:(forced_start_nc_time+5);
theor_slope_func = @(t) k/(1+sqrt(l))^2 * t; % k in min^{-1}, t in min
theor_slope_data = theor_slope_func(theor_time_mesh-forced_start_nc_time);
plot(theor_time_mesh, theor_slope_data, '-', 'LineWidth', 2, 'color', 'black');
%% Adding the maximal particle number
plot(xlim_nc_14, ones(1,2) * L/(l+sqrt(l)), ':', 'color', 'black',...
    'LineWidth', 2);
%% Adjusting the plot
xlim(xlim_nc_14);
ylim(ylim_nc_14);
xlabel('Time (min)');
ylabel('N_{pol}');
title(sprintf('Raw data, nc %i', nuc_cyc));


%% nc14 synced data
nuc_cyc=14;
% Loading data
input_filename = sprintf('%s_%s_nc_%i.mat', gene_name, dataset_name, nuc_cyc);
input_full_path = strcat(input_folder, input_filename);
load (input_full_path);
% Plotting
subplot(2, 2, 4);
hold on;
for i=1:ms2_spots_count
    if ms2_combined(i).nuc_cyc ~= nuc_cyc; continue; end;
    % Taking only those ms2 for which the calculated intersect was > 0
    if intersct_array(i)>0
        plot ((ms2_combined(i).Frame * mins_per_frame) + (forced_start_nc_time - intersct_array(i)), ...
            ms2_combined(i).Fluo/fluo_per_polymerase);        
    end;
end;
%% Adding the theoretical prediction
theor_time_mesh = forced_start_nc_time:0.05:(forced_start_nc_time+5);
theor_slope_func = @(t) k/(1+sqrt(l))^2 * t; % k in min^{-1}, t in min
theor_slope_data = theor_slope_func(theor_time_mesh-forced_start_nc_time);
plot(theor_time_mesh, theor_slope_data, '-', 'LineWidth', 2, 'color', 'black');
%% Adding the maximal particle number
plot(xlim_nc_14, ones(1,2) * L/(l+sqrt(l)), ':', 'color', 'black',...
    'LineWidth', 2);
%% Adjusting the plot
xlim(xlim_nc_14);
ylim(ylim_nc_14);
xlabel('Time (min)');
ylabel('N_{pol}');
title(sprintf('Synced data, nc %i', nuc_cyc));


%% Saving the figure
output_filename_png = sprintf('01_%s_%s_Raw_and_synced_data.png', gene_name, dataset_name);
utput_full_path_png = strcat(output_figures_folder, output_filename_png);
display('Saving figure...');
h_fig.PaperPositionMode = 'auto';
saveas(h_fig, utput_full_path_png, 'png');
display('Figure saved');







