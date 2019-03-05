clc;clear; close all
main_path=fileparts(mfilename('fullpath'));
addpath(genpath(main_path));
cd(main_path)
%% loading Data
% file_name='squat_data';
file_name='squat_data_PCA';
% file_name='words_data';
% file_name='cricket_data';

load(file_name)
if ~exist ('list_labels')
    disp('-- Please put labeling data into "list_labels"')
    disp('-- If you use labeling, please set all list_labels(i)=1')
end

%% calculate DTW of the data
if ~exist('alldistances')
    disp('--Please set distance matrix --> "alldistances"')
    disp('--or')
    disp('--To compute pairwise distances set each motion sequence "i" as "data_all{i}.data"')    
    alldistances=DTW_calc(data_all,'DTW',0,1);
end
D_path=fullfile(main_path,'data/DTW_result');
save(D_path, 'alldistances')

disp('--If you have your own gram-matrix, save it as "K" and in the following location')
K_path='';
% K_path=fullfile(main_path,'data/gram_mat');

%% ========= Important paramters
params.algorithm='LC-NNKSC'    %NNKSC: Unsupervised algorithm     LC-NNKSC: Supervised algorithm
params.T=4;  % limit on |x|_0  and |a|_0   --> increase it until you get the desired reconstruction error
params.list_labels=list_labels;
params.train_ratio=1;  % amount of data to be taken as trainin data, between (0 1]
params.D_path=D_path;
params.K_path=K_path;      % path of the kernel matrix if available
params.iter_end=100;   %last iteration of the algorithm
params.LC_betta=0.1;    % supervised weight for LC-NNKS algorithm

if strcmp(params.algorithm,'LC-NNKSC') && params.train_ratio==1
    params.train_ratio=0.75;    
end

%% ==== Optional or detailed parameters 
params.MP=max(list_labels)*params.T;  % dictionary size: default-->[class types]*T 
                                      % , change it if you have more variations inside the classes
params.T_a=params.T;   %limit on |a|_0 , 0--> no bound
params.a_gamma=0.1;    % a_{ij} similarity weight --> influences a_i to use similar training data from Kernel
params.class_a=1;    % limitting each a_i to use only one class of data
params.early_stop=1; %% stop the optimization when it converges
params.display=0;    % display the progress figures
detail_parameters;
%% None Negative Kernel Sparse Coding
results=NNKSC_func(params);
% NNKSC_func

%% Analyzing the results
fprintf('\n ==================== Analyzing the Dictionary =====\n');
[DS,D_intakes,D_intakes_ind]=D_interpretability(results);
fprintf('\n ==================== TSNE Visualization =======\n');
TSNE_visual(results,params)

return
%% Reconstructing new data (if available)
if ~exist('Kny')
    error('Please provide the kernel matrices for the new data (Knn, Kny)')
end
[X_new, err_x]=nnkomp_all(results.Adc,Kyy,Kny,Knn,params.T);

%% Accuracy of mapping for new data (if available)
[accuracy_new,conf_mat, i_missed]=LC_classification(X_new,results.Adc,results.H_tr,labels_new);