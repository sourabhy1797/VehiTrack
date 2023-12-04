%% This is mainstream file for obfuscating the indvidual node and node in trajectory using Lapcian Noise

%% Notations
% R = Real Location
% O = Obfuscation Location Vector
% post_vectors = cell containing Posterior vector for each node
% all_post_vec = All Posterior Vectors for a particular node
% obf_loc_posterior = Posterior Vector for Specific Obfuscated Location

% Bayesian formula  =>  posterior_prob = ( obf_prob * prior_distribution_R ) / denomi
% where 
% denomi = sum ( (obfuscated matrix vector for real location R) * pR )
% R = Real Location
% obf_nodes = Obfuscated Location
% obf_prob = Probability from obfuscation matrix for OI in RLV
% posterior_prob = posterior probability
% prior_distribution_R = prior distribution for each node


%% Main stream Code
addpath('./Functions/');

opts = detectImportOptions('./Dataset/smaller_nodes_rome.csv');
opts = setvartype(opts, 'osmid', 'int64');
opts = setvartype(opts, 'x', 'double');
opts = setvartype(opts, 'y', 'double');
df_nodes = readtable('./Dataset/smaller_nodes_rome.csv', opts);
df_nodes_big = readtable('./Dataset/nodes_rome.csv', opts);
df_edges = readtable('./Dataset/edges_rome.csv');

load('.\Dataset\rome_taxicab\traject_coord_approx.mat');
load('.\Dataset\rome_taxicab\trajectories_approx.mat'); 

load('reachable_mat_e50.mat')
% load('obf_mat_e150.mat')
load('C:\Users\sy0378.STUDENTS\Desktop\transformer_threat\HMM\data\trajectory_obfuscation\e100\obf_mat_new.mat')

% Extract the relevant columns (Node_ID, Latitude, Longitude)
nodes_data = table2array(df_nodes_big(:, 1:3));
node_lat_x = df_nodes.x;
node_lon_y = df_nodes.y;

%% New Nodes Dataset (Xinpeng's Dataset)
new_nodes_df = readtable('./Dataset/node.csv');
newColumnOrder = {'node', 'lat', 'lng'};
new_nodes_df = new_nodes_df(:,newColumnOrder);
new_nodes_df = renamevars(new_nodes_df, 'lat', 'y');
new_nodes_df = renamevars(new_nodes_df, 'lng', 'x');


%% Making both bigger and smaller data consistent just to avoid confusions among dataset
columnsToKeep = {'y', 'x'};
nodes_data_smaller = df_nodes(:, columnsToKeep);
large_coord = new_nodes_df(:, 2:3); 
smaller_coord = df_nodes(:,2:3);

% Use ismember to find the indices of matching coordinates
[~, indices] = ismember(large_coord, smaller_coord, 'rows');

% Create a new dataset containing the matching nodes
matchingNodes = new_nodes_df(indices > 0, :);
new_df_nodes = matchingNodes;

new_nodes_data = table2array(matchingNodes);
nodes_x = new_df_nodes.x;
nodes_y = new_df_nodes.y;
% scatter(nodes_x, nodes_y , 10, z_vector(1, :), "filled"); 
% hold on; 

new_node_ids = new_nodes_data(:,1);

%% New Trajectories from new Dataset (Xinpeng's Dataset)
trajectoryDataset = readtable('./Dataset/matching_result.csv');
trajectoryStrings = trajectoryDataset.Node_list;
newTrajectories = cell(length(trajectoryStrings), 1);

for i = 1 : numel(trajectoryStrings)
    trajectoryArray = str2num(strrep(trajectoryStrings{i}, ',', ' '));
    newTrajectories{i,1} = trajectoryArray;
end


%% Dropping of extra nodes from the trajectories
updated_trajectories = {};
for i = 1 : length(newTrajectories)
% for i = 1 : 5
    traj = newTrajectories{i,1};
    isMember = ismember(traj, new_node_ids);
    updated_traj = traj(isMember);
    if(~isempty(updated_traj))
        updated_trajectories{end+1, 1} = updated_traj;
    end
end

new_traj = {};
for i = 1:length(updated_trajectories)
    if(length(updated_trajectories{i})>1)           % Dropping all the trajectories which are of length 1
        new_traj{end+1,1} = updated_trajectories{i};
    end
end



%% Nodes Count in Trajectory
node_count = zeros(size(new_node_ids));     % Initialize a vector to store the counts
totalNodesEntireTrajectories = 0;

for i = 1:numel(new_traj)
    current_trajectory = new_traj{i};

    for j = 1:numel(current_trajectory)
        node_id = current_trajectory(j);
        totalNodesEntireTrajectories = totalNodesEntireTrajectories + 1;

        if ismember(node_id, new_node_ids)
            index = find(new_node_ids == node_id);
            node_count(index) = node_count(index) + 1;
        end
    end
end

prior_distribution = node_count/totalNodesEntireTrajectories;        % Stores Prior Distribution


%% Calculation of Obfuscated Node for each nodes based its occurence in trajectory

obfuscated_traj = {};
% i=30;
traj_set = {};
for a = 1:20
    obfuscated_traj = {};
    for i = 1 : 100
        disp("*****************")
        disp("Current Batch")
        disp(a)
        disp("Total Trajectories:")
        disp("100")
        disp("Current Trajectory:")
        disp(i)
        disp("------------------")
        traj = new_traj{i,1};
        obf_traj = {};
        for j = 1 : length(traj)
            id = traj(j);
            idx = find(new_node_ids == id);
            z_vec = z_vector(idx, :);
            obf_node = obfuscation1(new_df_nodes, z_vec);      % obfuscation of node
    %         obf_traj{end+1} = obf_node;
    %         obf_idx = find(new_node_ids == obf_node);
            obf_traj{end+1} = new_node_ids(obf_node);
    
        end
        obfuscated_traj{end+1,1} = obf_traj;
    end
    traj_set{end+1,1} = obfuscated_traj;
end

%% Calculation of Obfuscated node for each node in the dataset
% 
% obf_node_data = [];
% obf_x = [];
% obf_y = [];
% % i = 1;
% for i = 1 : length(new_node_ids)
%     disp(i)
%     disp(new_node_ids(i))
%     z_vec = z_vector(i, :);
%     obf_vec = z_vec;
%     obf_node = obfuscation100(new_df_nodes, z_vec);      % obfuscation of node
%     idx = new_node_ids(obf_node);
%     obf_node_data = [obf_node_data, idx];
%     obf_x = [obf_x, new_nodes_data(obf_node,3)];
%     obf_y = [obf_y, new_nodes_data(obf_node,2)];
% end
% obf_node_data = obf_node_data';

%% Ploting the Last trajectory

% obf_x = [];
% obf_y = [];
% real_x = [];
% real_y = [];
% for i = 1 : length(new_traj{100,1})
%     real_node_id = new_traj{100,1}(i);
%     idx = find(new_node_ids == real_node_id);
%     real_x = [real_x, new_df_nodes{idx, "x"}];
%     real_y = [real_y, new_df_nodes{idx, "y"}];
% 
%     obf_node_id = obfuscated_traj{100,1}(i);
%     obf_idx = find(new_node_ids == obf_node_id{1});
%     obf_x = [obf_x, new_df_nodes{obf_idx, "x"}];
%     obf_y = [obf_y, new_df_nodes{obf_idx, "y"}];
% 
% 
% end
% 
% scatter(obf_x, obf_y, 'o');
% hold on;
% scatter(real_x, real_y, '*')