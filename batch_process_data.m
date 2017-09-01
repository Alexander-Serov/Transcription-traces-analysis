%% This script launches the standard analysis on all nuclear cycles, genes and constructs


clear;

%% Constants
constants;



%% Initialize
set (0, 'DefaultAxesFontSize', 20);
gene_count = length(gene_names_array);
dataset_count = length(dataset_names_array);

for gene_ind = 1:gene_count
	for dataset_ind = 1:dataset_count
		for nuc_cyc = nuc_cyc_array
			analysis_7_new_data (gene_ind, dataset_ind, nuc_cyc);
		end;
	end;
end;


%% Plotting figures for the article
plot_mean_evolution_article;
plot_slopes_histogram_article;
plot_slopes_vs_AP_article_1_plot;
