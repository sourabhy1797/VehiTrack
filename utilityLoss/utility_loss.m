
%% This is mainstream file for generating the 1 posteriors for 1 node in trajectory

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

load('C:\Users\sy0378.STUDENTS\Dropbox\Sourabh Yadav\HMM\data\trajectory_obfuscation\e150\obf_mat.mat')
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
new_nodes_df.node(new_nodes_df.node == 0) = 50000;

new_edges_df = readtable('./Dataset/edge_weight.csv');
new_edges_df.s_node(new_edges_df.s_node == 0) = 50000;
new_edges_df.e_node(new_edges_df.e_node == 0) = 50000;



%% Graph Preparation

G = graph();
G = addnode(G, new_nodes_df);
numEdges = size(new_edges_df, 1);
for i = 1:numEdges
    s_node = new_edges_df.s_node(i);
    e_node = new_edges_df.e_node(i);
    max_speed = new_edges_df.max_speed(i);
    edge_length = new_edges_df.length(i);
    
    edge_weight = edge_length / max_speed;
    
    G = addedge(G, s_node, e_node, edge_length);
end


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

% obfuscated_traj = {};
% % i=30;
% for i = 1 : 100
%     disp("*****************")
%     disp("Total Trajectories:")
%     disp(length(new_traj))
%     disp("Current Trajectory:")
%     disp(i)
%     disp("------------------")
%     traj = new_traj{i,1};
%     obf_traj = {};
%     for j = 1 : length(traj)
%         id = traj(j);
%         idx = find(new_node_ids == id);
%         z_vec = z_vector(idx, :);
%         obf_node = obfuscation1(new_df_nodes, z_vec);      % obfuscation of node
% %         obf_idx = find(new_node_ids == id);
%         obf_traj{end+1} = new_node_ids(obf_node);
% 
%     end
%     obfuscated_traj{end+1,1} = obf_traj;
% end

% load('./Dataset/obfuscated_traj.mat')
load('C:\Users\sy0378.STUDENTS\Dropbox\Sourabh Yadav\HMM\data\trajectory_obfuscation\e50\obfuscated_traj.mat')

%% Utility Loss

pois = readtable('./Dataset/pois.csv');
pois = table2array(pois(:,1));
utility_loss_traj_whole = {};
new_pois = intersect(new_node_ids, pois);
new_pois_idx = find(ismember(new_node_ids, new_pois));
% For Nodes
obf_node_data = [];
obf_x = [];
obf_y = [];
% i = 1;
utility_loss_mat = z_vector;
% i =2;

for i = 1 : length(new_node_ids)
    real_id = new_node_ids(i);
    for j = 1 : length(new_node_ids)
        disp("i = ")
        disp(i)
        disp("j = ")
        disp(j)
        obf_id = new_node_ids(j);
        real_dist = [];
        obf_dist = [];
        for k = 1 : length(new_pois)
            [obf_path, obf_distance] = shortestpath(G, new_pois(k), obf_id);
            [real_path, real_distance] = shortestpath(G, new_pois(k), real_id);
            real_dist = [real_dist, real_distance];
            obf_dist = [obf_dist, obf_distance];
        end
        dist_diff = abs(obf_dist - real_dist);
        utility_loss_mat(i,j) = (nansum(dist_diff)/length(dist_diff));
    end
end

