



%% Constants
raw_data_plot_every_which_point = 1;
line_transparency = 1;


%% To make a legend, find now first curve from each bin
first_trace_in_bin_index = zeros(1, bins_count);
for bin = 1:bins_count
    
    index = find(ms2_index_to_bin_table == bin, 1, 'first');
    if index > 0  
        first_trace_in_bin_index(bin) = index;
    end;
end;

% Dropping zeros from the newly created variable
first_trace_in_bin_index(first_trace_in_bin_index==0) = [];

% Reordering indices to be able to add an automatic legend
reordered_indiced = setdiff(1:ms2_spots_count, first_trace_in_bin_index);
reordered_indiced = [first_trace_in_bin_index, reordered_indiced];



fig_hand = figure(12);
set_my_fig_size (fig_hand);
clf;
hold on;
for i=reordered_indiced
    if ms2_combined(i).nuc_cyc ~= nuc_cyc; continue; end;
    
    % Taking only those ms2 for which the calculated intersect was > 0
    if intersct_array(i)>0
        
        % Getting the current bin number
        current_bin = ms2_index_to_bin_table(i);
        if current_bin>bins_count || current_bin<=0
            continue;
        end;
        
        plot ((ms2_combined(i).Frame * mins_per_frame) + (forced_start_nc_time - intersct_array(i)), ...
            ms2_combined(i).Fluo/fluo_per_polymerase, 'color', [my_bins_color_sequence(current_bin, :), line_transparency],...
        'LineWidth', 1);
        
    end;

end;
% % alpha(0.9);
xlabel('Time, min');
ylabel('Number of active polymerases');
xlim(xlim_vector);
% xlim([47, 70]);
hold off;
% title(sprintf('Shifted fluorescence data for the %s construct', dataset_name));
% % 
% % % Adding the theoretical prediction
% % hold on;
% % theor_time_mesh = forced_start_nc_time:0.05:(forced_start_nc_time+4);
% % theor_slope_func = @(t) k/(1+sqrt(l))^2 * t; % k in min^{-1}, t in min
% % theor_slope_data = theor_slope_func(theor_time_mesh-forced_start_nc_time);
% % plot(theor_time_mesh, theor_slope_data, '--', 'LineWidth', 2, 'color', 'black');
% % 
% % 
% % % Adding the maximal particle number
% % temp_time_mesh = xlim_vector(1):0.05:xlim_vector(2);
% % plot(temp_time_mesh, ones(1,length(temp_time_mesh)) * L/(l+sqrt(l)), ':', 'color', 'black',...
% %     'LineWidth', 2);


% Adding a legend
legend_string = cell(bins_count, 1);
for i = 1:bins_count
    legend_string{i} = sprintf('AP = %.2f', bins_centers(i));
end;
legend(legend_string, 'Location', 'northwest');
hold off;



% Saving the figure
output_filename = sprintf('01b_%s_raw_shifted_data_%s_nc_%i.png', gene_name, dataset_name, nuc_cyc);
output_full_path = strcat(output_figures_folder, output_filename);
if bl_save_figures
    fig_hand.PaperPositionMode = 'auto';
    saveas(fig_hand, output_full_path, 'png');
    display('Figure saved');
end;








