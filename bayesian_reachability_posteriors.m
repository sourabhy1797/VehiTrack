%% This is mainstream file for generating the 1 posteriors for 1 node in trajectory

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

load('.\obf_matrix\obf_mat_e7_5.mat')
load('.\reachability_mat\reachable_mat.mat')
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
% 
% 
% 
% %% Nodes Count in Trajectory
% node_count = zeros(size(new_node_ids));     % Initialize a vector to store the counts
% totalNodesEntireTrajectories = 0;
% 
% for i = 1:numel(new_traj)
%     current_trajectory = new_traj{i};
% 
%     for j = 1:numel(current_trajectory)
%         node_id = current_trajectory(j);
%         totalNodesEntireTrajectories = totalNodesEntireTrajectories + 1;
% 
%         if ismember(node_id, new_node_ids)
%             index = find(new_node_ids == node_id);
%             node_count(index) = node_count(index) + 1;
%         end
%     end
% end
% 
% prior_distribution = node_count/totalNodesEntireTrajectories;        % Stores Prior Distribution
% 
% 
% %% Calculation of Posteriors for each nodes based its occurence in trajectory
% 
% % for k1 = 1 : 20
% post_vectors = {};
% for i = 1 : 100
%     % for i = 1 : length(new_traj)
%     disp("*****************")
%     disp("Total Trajectories:")
%     disp(length(new_traj))
%     disp("Current Trajectory:")
%     disp(i)
%     disp("------------------")
%     traj = new_traj{i,1};
%     traj_posts = {};
%     loc_count = [];
%     for j = 1 : length(traj)
%         id = traj(j);
%         idx = find(new_node_ids == id);
%         z_vec = z_vector(idx, :);
%         obf_node = obfuscation_node(new_df_nodes, z_vec);      % obfuscation of node
%         obf_loc_posterior = [];
%         reach_vec = reachable_vec(:, obf_node);
%         Real_idx = find(reach_vec == 1);
%         denomi = [];
%         for k = 1:length(Real_idx)
%             R = Real_idx(k);                    % R = Real Location
%             obf_prob = z_vector(R, obf_node);                        % obf_prob = Probability from obfuscation matrix for OI in RLV
%             priorDistribution_R = prior_distribution(R);                    % prior_distribution_R = prior distribution for k_th node in R
%             denomi = [denomi, (obf_prob * priorDistribution_R)];                                  % Bayes formula's denominator
%         end
%         for l = 1:length(reachable_vec)
%             if reach_vec(l, 1) > 0
%                 obf_prob = z_vector(l, obf_node);
%                 priorDistribution_R = prior_distribution(l);
%                 posterior_prob = (obf_prob * priorDistribution_R)/sum(denomi);                                   % Bayes Formula
%                 obf_loc_posterior(1, l) = posterior_prob;
%             else
%                 obf_loc_posterior(1, l) = 0;
%             end
%         end
%         traj_posts{end+1} = obf_loc_posterior;
%     end
%     post_vectors{end+1, 1} = traj_posts;
% end
% % filename = './test/pv';
% % filename = strcat(filename, num2str(k1));
% % filename = strcat(filename, '.mat');
% % save(filename, 'post_vectors')
% % end
% 
% % filepath = 'C:\Users\sy0378.STUDENTS\Desktop\Updated Dataset\pv3\pv.mat';
% % save(filepath, 'prior_distribution');


