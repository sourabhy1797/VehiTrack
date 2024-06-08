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

% load('obf_mat_e150.mat')
load('.\obf_matrix\obf_mat_e5.mat')
load('.\reachability_mat\reachable_mat.mat')
load('.\Dataset\distanceMatrix.mat')

% Extract the relevant columns (Node_ID, Latitude, Longitude)
nodes_data = table2array(df_nodes_big(:, 1:3));
node_lat_x = df_nodes.x;
node_lon_y = df_nodes.y;

%% New Nodes Dataset (Xinpeng's Dataset)        (Changing the Node ID 0 to 50000)
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

% G = graph();
% G = addnode(G, new_nodes_df);
% numEdges = size(new_edges_df, 1);
% for i = 1:numEdges
%     s_node = new_edges_df.s_node(i);
%     e_node = new_edges_df.e_node(i);
%     max_speed = new_edges_df.max_speed(i);
%     edge_length = new_edges_df.length(i);
% 
%     edge_weight = edge_length / max_speed;
% 
%     G = addedge(G, s_node, e_node, edge_length);
% end



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


%% Computation of Violations
load('./Dataset/location_set_new100.mat')
new_loc_set = {};
for i = 1:length(location_set)
    traject = location_set{i,1};
    for j = 1:length(traject(:,1))
        loc_set = traject(j,:);
        isMember = ismember(loc_set, new_node_ids);
        updated_loc_set = loc_set(isMember);
        if(~isempty(updated_loc_set))
            new_loc_set{i,1}(j, :) = updated_loc_set;
        end
    end
end

obfuscation_matrix = z_vector;  % Your 1000x1000 obfuscation matrix
trajectories = new_traj;  % Trajectory of node IDs
% target_region_nodes = new_loc_set{1};  % Cell array containing target region nodes for each node in the trajectory
epsilon = 5;  % Privacy parameter
distances = distanceMatrix;  % Matrix containing distances between nodes



% num_trajectories = length(trajectories);
num_trajectories = 100;
num_nodes = size(obfuscation_matrix, 1);

average_percentage_per_node = zeros(num_nodes, 1);
total_violation_percentage = 0;

for t = 1:num_trajectories
    trajectory = trajectories{t};
    target_region_nodes = new_loc_set{t};
    % violation_percentages = compute_violation_percentages(obfuscation_matrix, trajectory, target_region_nodes{t}, epsilon, distances);
    num_nodes = size(obfuscation_matrix, 1);
    num_trajectory_nodes = length(trajectory);
    violation_percentages = zeros(num_trajectory_nodes, 1);

    for i = 1:num_trajectory_nodes
        current_node = trajectory(i);
        current_node_idx = find(new_node_ids == current_node);
        current_target_region_nodes = target_region_nodes(i,1:10);
        % Initialize an array to store the indexes
        current_target_region_nodes_indexes = zeros(1, numel(current_target_region_nodes));

        % Loop through each node in the sublist
        for j = 1:numel(current_target_region_nodes)
            % Find the index of the current node in the main list
            current_target_region_nodes_indexes(j) = find(new_node_ids == current_target_region_nodes(j));
        end
        num_current_target_nodes = length(current_target_region_nodes_indexes);
        num_non_target_nodes = num_nodes - num_current_target_nodes;

        non_target_nodes = setdiff(1:num_nodes, current_target_region_nodes_indexes);

        non_target_probability_sum = sum(obfuscation_matrix(current_node_idx, non_target_nodes));
        target_probability_sum = sum(obfuscation_matrix(current_node_idx, current_target_region_nodes_indexes));

        % Calculate violation for each target region node
        for j = 1:num_current_target_nodes
            target_node = current_target_region_nodes_indexes(j);
            distance = distances(current_node_idx, target_node);  % Distance between current node and target node

            if target_probability_sum > 0
                expected_probability = target_probability_sum * exp(distance * epsilon);
                if non_target_probability_sum > expected_probability
                    violation_percentages(i) = violation_percentages(i) + (non_target_probability_sum - expected_probability) / num_non_target_nodes;
                end
            end
        end

        violation_percentages(i) = violation_percentages(i) * 100;  % Convert to percentage
    end
    
    % Accumulate violation percentages for each node
    for n = 1:length(trajectory)
        current_node_idx = find(new_node_ids == trajectory(n));
        average_percentage_per_node(current_node_idx) = average_percentage_per_node(current_node_idx) + violation_percentages(n);
    end
    
    % Accumulate total violation percentage
    total_violation_percentage = total_violation_percentage + sum(violation_percentages);
end

% Compute average violation percentage for each node
average_percentage_per_node = average_percentage_per_node / num_trajectories;

% Compute average violation percentage for all trajectories combined
average_percentage_all_trajectories = total_violation_percentage / sum(cellfun(@length, trajectories));

disp('Average violation percentage for each node:');
disp(average_percentage_per_node);

disp('Average violation percentage for all trajectories combined:');
disp(average_percentage_all_trajectories);