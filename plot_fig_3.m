

%% Plotting 2D curves in 3D. Np_mean as a function of time and AP position






% bin_fluo_data_shifted_mean (frame, bin);


% Choose time frame
% frame = 52*60;
% Choose bin
fill_opacity = 0.3;


%% Initialization
% Finding bins that are not empty
bin_not_empty = sum(~isnan(bin_fluo_data_shifted_mean),1)>0;
not_empty_bin_indices = find(bin_not_empty>0);





%% Plotting a 2D slice in 3D
fig_hand = figure(13);
set_my_fig_size (fig_hand);
clf;
hold on;

%% Plotting the lines
for bin=not_empty_bin_indices
    AP_constant_mesh = ones(1, length(new_time_mesh)) * bins_centers(bin);
    
    plot3(new_time_mesh, AP_constant_mesh, bin_fluo_data_shifted_mean(:, bin)/fluo_per_polymerase,...
        'color', my_bins_color_sequence(bin, :), 'LineWidth', 2);
    
    
end;

%% Plotting the fill
for bin=not_empty_bin_indices
    
    %% Preparing vertices to make the 3D error region fill
    t_vertices = new_time_mesh(1:bootstrap_only_each_frame:end);
    t_vertices = [t_vertices, flip(t_vertices)];
    
    AP_vertices = ones(1, length(t_vertices)) * bins_centers(bin);
    
    fluo_vertices = [...
        bin_fluo_data_shifted_mean(1:bootstrap_only_each_frame:end, bin) + bin_fluo_data_shifted_confidence_intervals(1:bootstrap_only_each_frame:end, bin,1);...
        flip(bin_fluo_data_shifted_mean(1:bootstrap_only_each_frame:end, bin) - bin_fluo_data_shifted_confidence_intervals(1:bootstrap_only_each_frame:end, bin,2))]/fluo_per_polymerase;
    
    %% Some of the fluo vertices can be NaN. So I need to remove those indices from all three vectors
    not_nan_indices = ~isnan(fluo_vertices);
    t_vertices = t_vertices(not_nan_indices);
    AP_vertices = AP_vertices(not_nan_indices);
    fluo_vertices = fluo_vertices(not_nan_indices);
    
    
    %% Creating the fill
    fill3(t_vertices, AP_vertices, fluo_vertices', my_bins_color_sequence(bin, :), 'EdgeColor', 'None',...
        'FaceAlpha', fill_opacity);
    
   
    
end;

% Adding a legend
legend_string = cell(sum(bin_not_empty), 1);
counter = 1;
for i = not_empty_bin_indices
    legend_string{counter} = sprintf('AP = %.2f', bins_centers(i));
    counter = counter + 1;
end;
legend(legend_string, 'Location', 'northeast');


hold off;
    
xlabel('Time (min)');
ylabel('AP position');
zlabel('$\overline{N_p}$', 'interpreter', 'latex');

y_lim_vector = [bins_borders(1), bins_borders(end)];
xlim(xlim_vector);
ylim(y_lim_vector);

% title(sprintf('Time and AP position dependence for %s construct', dataset_name));

view(20,40);


% Saving the figure
output_filename = sprintf('03_%s_3D_time_AP_dependence_%s_nc_%i.png', gene_name, dataset_name, nuc_cyc);
output_full_path = strcat(output_figures_folder, output_filename);
if bl_save_figures
    fig_hand.PaperPositionMode = 'auto';
    saveas(fig_hand, output_full_path, 'png');
    display('Figure saved');
end;







