%% Load and pre-sync fluorescence data from different embryos for the given gene, construct and nuclear cycle


function [Data, datasets, dataset_count, time_alignment_mins] = choose_datasets_and_time_alignment(gene_name, dataset_name, nuc_cyc)

%% Constants
constants;



%% Initialize
Data = [];



%% Load and align data
% The time alignment constants are stored directly here

if nuc_cyc == 14
    %%
    if strcmp(gene_name, 'HunchBack') && strcmp(dataset_name, 'bac')
        load(strcat(data_folder, gene_name, '_', dataset_name));
        dataset_count = length(Data);    
        datasets = 1:11;    % 1:11
        % Do not use: 
            time_alignment_mins = [ 29, 18.5, 12.5, 0.5, 13, ...
                14, 7.5, 14, 23.5, 23, ...
                2.5];
                                
% % %             time_alignment_mins = [ 29.5, 19, 12.5, 5.5, 16.5, ...
% % %                 18.5, 12.5, 18, 27, 26, ...
% % %                 6];

    elseif strcmp(gene_name, 'HunchBack') && strcmp(dataset_name, 'no_primary')
        load(strcat(data_folder, gene_name, '_', dataset_name));
        dataset_count = length(Data);    
        datasets = [1:8, 11];
    % %     % Do not use: 9, 10, 12, 13, 14, 15, 16, 17
        time_alignment_mins = [ 32.5, 30, 16.5, 11, 24.5, ...
                                    18, 10.5, 10.5, 24.5, 11, ...
                                    0.5, 24.5, 25.5, 26.5, 26.5, ...
                                    0, 24.5];   


    elseif strcmp(gene_name, 'HunchBack') && strcmp(dataset_name, 'no_shadow')
        load(strcat(data_folder, gene_name, '_', dataset_name));
        dataset_count = length(Data);    
        datasets = 1:8;
    % %     % Do not use: 9, 10, 11, 12
            time_alignment_mins = [ 0, 6, 27, 16, 5.5, ...
                                    10, 27.5, 22, 0, 0, ...
                                    0];                              

    elseif strcmp(gene_name, 'Knirps') && strcmp(dataset_name, 'bac')
        load(strcat(data_folder, gene_name, '_', dataset_name));
        dataset_count = length(Data);    
        datasets = 1:5;
    % %     % Do not use:
            time_alignment_mins = [ 9.5, 7.5, -8, 3, 0];    

    elseif strcmp(gene_name, 'Knirps') && strcmp(dataset_name, 'no_primary')
        load(strcat(data_folder, gene_name, '_', dataset_name));
        dataset_count = length(Data);    
        datasets = 1:2;
    % %     % Do not use: 3
            time_alignment_mins = [ -1, -0.5, 0];                               

    elseif strcmp(gene_name, 'Knirps') && strcmp(dataset_name, 'no_shadow')
        load(strcat(data_folder, gene_name, '_', dataset_name));
        dataset_count = length(Data);    
        datasets = 4;
    % %     Unreliable data sets: 1, 2 
            time_alignment_mins = [ 22.5, 32.5, 35.5, 21];                      

    elseif strcmp(gene_name, 'SNAIL') && strcmp(dataset_name, 'bac')
        load(strcat(data_folder, gene_name, '_', dataset_name));
        dataset_count = length(Data);    
        datasets = 1:5;
    % %     % Do not use: 
            time_alignment_mins = [ 6, 45, 20, 34.5, 34];      

    elseif strcmp(gene_name, 'SNAIL') && strcmp(dataset_name, 'no_primary')
        load(strcat(data_folder, gene_name, '_', dataset_name));
        dataset_count = length(Data);    
        datasets = 1:4;
    % %     % Do not use: 
            time_alignment_mins = [ 13.5, 38.5, 26.9, 13];    

    elseif strcmp(gene_name, 'SNAIL') && strcmp(dataset_name, 'no_shadow')
        load(strcat(data_folder, gene_name, '_', dataset_name));
        dataset_count = length(Data);    
        datasets = 1;
    % %     % Do not use: 
            time_alignment_mins = [ 24];                                 


    end;
    
    
    
elseif nuc_cyc == 13
    %%
    if strcmp(gene_name, 'HunchBack') && strcmp(dataset_name, 'bac')
        load(strcat(data_folder, gene_name, '_', dataset_name));
        dataset_count = length(Data);    
        datasets = 1:11;
    % %     % Do not use: 
            time_alignment_mins = [ 29.5, 20, 14.5, 5.5, 16.5, ...
                                    18.5, 13, 18, 27, 26, ...
                                    6];

    elseif strcmp(gene_name, 'HunchBack') && strcmp(dataset_name, 'no_primary')
        load(strcat(data_folder, gene_name, '_', dataset_name));
        dataset_count = length(Data);    
        datasets = [1:15, 17];
    % %     % Do not use: 16
            time_alignment_mins = [ 31, 34, 17, 17, 26, ...
                                    20, 14.5, 14.5, 25.5, 12.5, ...
                                    2, 26, 25.5, 28, 26.5, ...
                                    0, 25];   


    elseif strcmp(gene_name, 'HunchBack') && strcmp(dataset_name, 'no_shadow')
        load(strcat(data_folder, gene_name, '_', dataset_name));
        dataset_count = length(Data);    
        datasets = 1:12;
    % %     % Do not use: 
            time_alignment_mins = [ 6, 10, 27.5, 19, 8.5, ...
                                    14, 28, 25.5, 14, 27, ...
                                    25, 27.5];                              

    elseif strcmp(gene_name, 'Knirps') && strcmp(dataset_name, 'bac')
        load(strcat(data_folder, gene_name, '_', dataset_name));
        dataset_count = length(Data);    
        datasets = 1:4;
    % %     % Do not use: 5
            time_alignment_mins = [ 13, 11, 5, 8, 0];    

            
    elseif strcmp(gene_name, 'Knirps') && strcmp(dataset_name, 'no_primary')
        load(strcat(data_folder, gene_name, '_', dataset_name));
        dataset_count = length(Data);    
        datasets = [];
    % %     % Do not use: 1,2, 3. There is just no data on no_primary in nc13
            time_alignment_mins = [ -1, -0.5, 0];                               

    elseif strcmp(gene_name, 'Knirps') && strcmp(dataset_name, 'no_shadow')
        load(strcat(data_folder, gene_name, '_', dataset_name));
        dataset_count = length(Data);    
        datasets = [];
    % %     Do not use: 1 -- 5. Just no data.
            time_alignment_mins = [ 22.5, 32.5, 35.5, 21];                      

    elseif strcmp(gene_name, 'SNAIL') && strcmp(dataset_name, 'bac')
        load(strcat(data_folder, gene_name, '_', dataset_name));
        dataset_count = length(Data);    
        datasets = [1, 3:5];
    % %     % Do not use: 2 (no data), 
            time_alignment_mins = [ 13, 25, 26, 34.5, 34];      

    elseif strcmp(gene_name, 'SNAIL') && strcmp(dataset_name, 'no_primary')
        load(strcat(data_folder, gene_name, '_', dataset_name));
        dataset_count = length(Data);    
        datasets = [1, 3, 4];
    % %     % Do not use: 2 (no data), 
            time_alignment_mins = [ 23.5, 38.5, 32, 16];    

    elseif strcmp(gene_name, 'SNAIL') && strcmp(dataset_name, 'no_shadow')
        load(strcat(data_folder, gene_name, '_', dataset_name));
        dataset_count = length(Data);    
        datasets = 1;
    % %     % Do not use: --
            time_alignment_mins = [ 27.5];                                 


    end;
    
end;










