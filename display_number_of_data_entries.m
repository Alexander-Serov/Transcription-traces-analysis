

%% Constants
dataset_names_array = {'bac', 'no_primary', 'no_shadow'};
dataset_short_names_array = {'bac', 'no_pr', 'no_sh'};
gene_names_array = {'HunchBack', 'Knirps', 'SNAIL'};
gene_short_names_array = {'Hb', 'Kn', 'Sn'};
nuc_cyc_array = [13, 14];
input_folder = './processed_data/';


%% Initialization
fprintf('\n');


for gene_index = 1: length(gene_names_array)
    gene_name = gene_names_array{gene_index};
    for construct_index = 1:length(dataset_names_array)
        dataset_name = dataset_names_array{construct_index};    
        number_of_not_nan_points = zeros(1,2);
        for nuc_cyc_index = 1:length(nuc_cyc_array)
            nuc_cyc = nuc_cyc_array(nuc_cyc_index);
            % Loading data
            input_filename = sprintf('%s_%s_nc_%i.mat', gene_name, dataset_name, nuc_cyc);
            input_full_path = strcat(input_folder, input_filename);
            load (input_full_path);
            number_of_not_nan_points(nuc_cyc_index) = sum(~isnan(normalized_slopes_array));

        end;
        fprintf('Number of analyzed traces for %s (%s):\tnc13 - %i\tnc14 - %i\ttotal - %i\n',...
            gene_name, dataset_name, number_of_not_nan_points(1), number_of_not_nan_points(2),...
            sum(number_of_not_nan_points));
    end;
end;









