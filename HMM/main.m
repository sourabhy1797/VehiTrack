addpath('./func'); 
addpath('./func/haversine/'); 


%% HMM model
% load('.\data\obf_matrix.mat'); 
% load('.\data\obfuscated_traj.mat'); 
% load('.\data\new_traj.mat'); 
% new_nodes_df = csvread('.\data\new_nodes_df.csv'); 

% load('.\data\trajectory_obfuscation\e100\new_obf_nodes.mat');
load('.\data\trajectory_obfuscation\e50\obf_mat.mat'); 
load('.\data\trajectory_obfuscation\e50\obfuscated_traj.mat'); 
load('.\data\trajectory_obfuscation\e50\new_traj.mat'); 
new_nodes_df = csvread('.\data\new_nodes_df.csv'); 
matrix = readmatrix('.\data\transition.xlsx'); 
load('.\data\trans_matrix.mat'); 

node_id = new_nodes_df(:, 1)+1;                                             % The node ID starts from 0, while in the transition matrix, the index start from 1. 
                                                                            % To retrieve the corresponding submatrix, we need this index conversion

trans_matrix_ = trans_matrix(node_id, node_id); 
trans_matrix_ = full(trans_matrix_); 

EIE_HMM = zeros(100, 1); 

error_ = zeros(100, 1); 

trans_matrix_ = trans_matrix_ + 0.000000001; 

for i = 1:1:size(trans_matrix_, 1)
    trans_matrix_(i, :) = trans_matrix_(i, :)/sum(trans_matrix_(i, :)); 
end

states = 1:size(node_id, 1); 


%% HMM 
for i = 1:1:100
    obf_loc = []; 
    for j = 1:1:size(obfuscated_traj{i}, 2)
        obf_loc = [obf_loc, find(new_nodes_df == obfuscated_traj{i}{j})]; 
    end
    real_loc = []; 
    for j = 1:1:size(new_traj{i}, 2)
        real_loc = [real_loc, find(new_nodes_df == new_traj{i}(j))]; 
    end

    est_loc = hmmviterbi(obf_loc, trans_matrix_, z_vector);

    for j = 1:1:size(obfuscated_traj{i}, 2)
        real_x(1, j) = new_nodes_df(real_loc(1, j), 2);
        real_y(1, j) = new_nodes_df(real_loc(1, j), 3);
        obf_x(1, j) = new_nodes_df(obf_loc(1, j), 2);
        obf_y(1, j) = new_nodes_df(obf_loc(1, j), 3);


        est_x(1, j) = new_nodes_df(est_loc(1, j), 2);
        est_y(1, j) = new_nodes_df(est_loc(1, j), 3);

        
        [km, ~,  ~] = haversine([est_x(1, j), est_y(1, j)], [real_x(1, j), real_y(1, j)]); 
        EIE_HMM(i, 1) = EIE_HMM(i, 1) + km; 

        % Testing
        [km, ~,  ~] = haversine([obf_x(1, j), obf_y(1, j)], [real_x(1, j), real_y(1, j)]);

        error_(i, 1) = error_(i, 1) + km;

    end
    EIE_HMM(i, 1) = EIE_HMM(i, 1)/size(obfuscated_traj{i}, 2); 
    error_(i, 1) = error_(i, 1)/size(obfuscated_traj{i}, 2);
end
