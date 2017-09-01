%% This function stores plotting functions from the analysis code
% I think this file is not currently used


%%%%% ===== PLOT 1 =====

%% Plotting the evolution of the mean fluorescence
figure(3);
clf;
subplot(1,2,1);
hold on;
subplot(1,2,2);
hold on;
str_legends=cell(1, bins_count);

for bin=1:bins_count
    % Non-shifted data
    subplot(1,2,1);
    plot(original_time_mesh,...
        bin_fluo_data_non_shifted_mean(:,bin)/fluo_per_polymerase, 'LineWidth', 2);
    
    % Shifted data
    subplot(1,2,2);
    if bin==bootstrap_only_bin_number
        
        error_handle = errorbar(time_mesh(1:bootstrap_only_each_frame:end),...
            bin_fluo_data_shifted_mean(1:bootstrap_only_each_frame:end,bin)/fluo_per_polymerase,...
            bin_fluo_data_shifted_confidence_intervals(1:bootstrap_only_each_frame:end, bin,1)/fluo_per_polymerase,...
            bin_fluo_data_shifted_confidence_intervals(1:bootstrap_only_each_frame:end, bin,2)/fluo_per_polymerase,...
        '-o', 'LineWidth', 2);
%         errorbar_tick(error_handle, 50);
    else
        plot(time_mesh,...
        bin_fluo_data_shifted_mean(:,bin)/fluo_per_polymerase, 'LineWidth', 2);
    end;
    
    % Saving the legend
    str_legends{bin} = sprintf('AP = %.2f', bins_borders(bin) + bin_width/2);
    
end;

% xlim_vector = [48, 60]; %[48, 56]
ylim_vector = [0,110];

subplot(1,2,1);
xlabel('Time, min');
ylabel('Number of polymerases');
xlim(xlim_vector);
ylim(ylim_vector);
title('Non-shifted data');
legend(str_legends, 'Location', 'northwest');

% Adding the fill
x_vertices = forced_start_nc_time + [0, 20, 20, 0];
y_vertices = [0,...
    k/(1+sqrt(l_min))^2 * (x_vertices(2) - x_vertices(1)),...
    k/(1+sqrt(l_max))^2 * (x_vertices(3) - x_vertices(1)),...
    0];


subplot(1,2,2);
xlabel('Time, min');
ylabel('Number of polymerases');
xlim(xlim_vector);
ylim(ylim_vector);
title('Shifted data');
legend(str_legends, 'Location', 'northwest');
fill(x_vertices, y_vertices, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'None');

% Adding the theoretical prediction to the second plot
subplot(1,2,2);
theor_time_mesh = forced_start_nc_time + (0:time_step:init_slope_length);
theor_slope_func = @(t) k/(1+sqrt(l))^2 * t; % k in min^{-1}, t in min
theor_slope_data = theor_slope_func(theor_time_mesh-forced_start_nc_time);

plot(theor_time_mesh, theor_slope_data, '--', 'LineWidth', 2, 'color', 'black');

% Adding the maximal particle number
temp_time_mesh = xlim_vector(1):time_step:xlim_vector(2);
plot(temp_time_mesh, ones(1,length(temp_time_mesh)) * L/(l+sqrt(l)), ':', 'color', 'black',...
    'LineWidth', 2);


%% Plotting variance as a function of mean as suggested by Tkacik
figure(7);
clf;

for bin=1:bins_count
    % Non-shifted data
    subplot(1,2,1);
    scatter(bin_fluo_data_non_shifted_mean(:,bin)/fluo_per_polymerase,...
        (bin_fluo_data_non_shifted_STD(:,bin)/fluo_per_polymerase).^2, 'LineWidth', 2);
    hold on;
    
    % Shifted data
    subplot(1,2,2);
    scatter(bin_fluo_data_shifted_mean(:,bin)/fluo_per_polymerase,...
        (bin_fluo_data_shifted_STD(:,bin)/fluo_per_polymerase).^2, 'LineWidth', 2);
    hold on;
end;

ylim_vector = [0, 500];

subplot(1,2,1);
xlabel('<N>');
ylabel('Var(N)');
% xlim(xlim_vector);
ylim(ylim_vector);
title('Non-shifted data');
hold off;

subplot(1,2,2);
xlabel('<N>');
ylabel('Var(N)');
% xlim(xlim_vector);
ylim(ylim_vector);
title('Shifted data');
hold off;



%% Plotting the number of spots in a bin
figure(5);
clf;
subplot(1,2,1);
hold on;
subplot(1,2,2);
hold on;

for bin=1:bins_count
    % Non-shifted data
    subplot(1,2,1);
    plot(original_time_mesh,...
        bin_fluo_data_non_shifted_count(:,bin), 'LineWidth', 2);
    
    % Shifted data
    subplot(1,2,2);
    plot(time_mesh,...
        bin_fluo_data_shifted_count(:,bin), 'LineWidth', 2);
end;

% xlim_vector = [48, 60]; % 56
% % ylim_vector = [0,25];

subplot(1,2,1);
xlabel('Time, min');
ylabel('Number of data points');
xlim(xlim_vector);
% ylim(ylim_vector);
title('Non-shifted data');

subplot(1,2,2);
xlabel('Time, min');
ylabel('Number of data points');
xlim(xlim_vector);
% ylim(ylim_vector);
title('Shifted data');





% 
% %% Looking how these fluctuations are distributed across the ensemble
% figure(6);
% clf;
% % % subplot(1,2,1);
% % % hold on;
% % % subplot(1,2,2);
% % % hold on;
% 
% selected_new_mesh_frame_index = 847;
% % hold on;
% for bin=1:bins_count
%     subplot(1, bins_count, bin);
%     histogram(bin_fluo_data_shifted{selected_new_mesh_frame_index, bin}/fluo_per_polymerase, 1:8:80,...
%         'Normalization', 'count');
%     title(sprintf('AP = %.2f. T = %.0f min', ...
%         bins_borders(bin) + bin_width/2,...
%         new_time_mesh(selected_new_mesh_frame_index)));
%     xlabel('N');
%     ylabel('PDF');
% end;



%% Plotting data distribution in bins
% figure(8);




%% Now I want to check how changes with time the number of data I calculate
% means and stds over


