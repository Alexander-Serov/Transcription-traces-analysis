


%% Constants

secs_per_frame = 37;  % in seconds
mins_per_frame = secs_per_frame/60;  % in seconds
fluo_per_polymerase = 3.5e4/100; % in Stucken
forced_start_nc_time_13 = 35;
forced_start_nc_time_14 = 50;
manual_bins=true;
nuc_cyc = 14;
time_shift_limit = 5;  % in mins
k = 26*60;     % Bulk jump rate in bp/min
dk = 2 * 60;   % Error in k, bp/min
L = 5400;      % Gene length in base pairs
l = 45;        % Polymerase footprint in bp
dl = 5;         % Error in l
l_min = 40;     % Min polymerase footprint in bp
l_max = 60;     % Max polymerase footprint in bp
tau_exp = 5;    % in Seconds, for theoretical predictions
L_array = [5400, 1668, 3049];       % Gene lengths in base pairs for Hb, Sn and Kn.
dL_array = [0, 0, 0];         % Error in gene length for Hb, Sn and Kn.
bootstrap_samples_count = 2000;
bootstrap_only_bin_value = 0.37;        % Only used for plots. Bootstrap calculated for all bins
bootstrap_only_each_frame = 5;  % 5, 8, 12
% N_filter_threshold = 15;
integral_threshold_value = 200;     % in pol * mins through the nc14. Filters out traces that do not get above a certain value.
% filter_time_interval_mins = [48, 60];
init_slope_length = 2.5;      % in mins
raw_data_plot_every_which_point = 100;
output_figures_folder = './figures_for_article/';
bl_save_figures = false;
% bl_stop_early = true;
bl_stop_early = false;
half_width_max_rgn_ind = 2; % Number of points to consider around maximum when calculating the value of
                            % the plateau of concentration
input_folder = './processed_data/';
dataset_name_array = {'bac', 'no_primary', 'no_shadow'};
gene_names_array = {'HunchBack', 'SNAIL', 'Knirps'};
nuc_cyc_array = [13, 14];