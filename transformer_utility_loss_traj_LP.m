
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

% load('obf_mat_e150.mat')
load('.\obf_matrix\obf_mat_e10.mat')
load('.\reachability_mat\reachable_mat.mat')

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
pois = readtable('./Dataset/pois.csv');
pois = table2array(pois(:,1));
new_pois = intersect(new_node_ids, pois);
new_pois_idx = find(ismember(new_node_ids, new_pois));

load("./Dataset/ranges_nodes.mat")
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

obfuscated_traj = {};
obf_mat_traj = {};
utility_loss_mat = {};
expected_utility_loss_alltrajs = {};
average_expected_utility_loss = [];
post_vectors = {};
% i=13;
count = 0;
for i = 1 : 100
    disp("Current Trajectory:")
    disp(i)
    disp("------------------")
    traj = new_traj{i,1};
    obf_traj = {};
    node_utility_loss_mat = {};
    node_obf_mat = {};
    expected_utility_loss_traj = [];
    traj_posts = {};
    prob_loc_set = new_loc_set{i};
    
    % j =19;
    for j = 1 : length(traj)
        id = traj(j);
        idx = find(new_node_ids == id);
        probable_loc = prob_loc_set(j,1:5);
%         nodes_in_range = range_nodes{idx, :};
        Real_idx = find(ismember(new_node_ids, probable_loc));
        nodes_in_range = new_node_ids(Real_idx);
%         reach_vec = reachable_vec(:, idx);
%         Real_idx = find(reach_vec == 1);
%         disp("Reachability Matrix Received")
% %         nodes_in_range = range_nodes{idx, :};
%         if numel(Real_idx) > 40
%             % Generate a random permutation of indices
%             randomIndices = randperm(numel(Real_idx));
% 
%             % Select the first 40 elements
%             selectedNodeIDs = Real_idx(randomIndices(1:30));
%         else
%             % If length is less than or equal to 40, keep the list same
%             selectedNodeIDs = Real_idx;
%         end
%         Real_idx = selectedNodeIDs;
%         if ~any(Real_idx == idx)
%             % If not present, add the specific number to the array
%             Real_idx = [Real_idx' idx]';
%         end
%         nodes_in_range = new_node_ids(Real_idx);
        utility_loss = zeros(length(nodes_in_range), length(nodes_in_range));
        
        % utility loss computation
        for m = 1 : length(nodes_in_range)
            real_id = nodes_in_range(m);
            for n = 1 : length(nodes_in_range)
                obf_id = nodes_in_range(n);
                z_vec = z_vector(n,:);
                z_vec_probs = z_vec(new_pois_idx);
                real_dist = [];
                obf_dist = [];
                for k = 1 : length(new_pois)
                    [obf_path, obf_distance] = shortestpath(G, new_pois(k), obf_id);
                    [real_path, real_distance] = shortestpath(G, new_pois(k), real_id);
                    real_dist = [real_dist, real_distance];
                    obf_dist = [obf_dist, obf_distance];
                end
                dist_diff = abs(obf_dist - real_dist);
                utility_loss(m,n) = (nansum(dist_diff));
            end
        end

        distanceMatrix = zeros(length(nodes_in_range), length(nodes_in_range));

        for o = 1 : length(nodes_in_range)
            for p = 1 : length(nodes_in_range)
                distanceMatrix(o,p) = haversine_calc(matchingNodes.y(o), matchingNodes.x(o), matchingNodes.y(p), matchingNodes.x(p));
            end
        end

        node_utility_loss_mat{end+1} = utility_loss;

        try
            obf_mat = LPObfuscationMatrx_test(utility_loss, distanceMatrix, 5);
            node_obf_mat{end+1} = obf_mat;
            node_ids_in_range_mat = find(nodes_in_range == id);
            z_vec = obf_mat(node_ids_in_range_mat, :);
        catch
            % loc_set_z_vec = zeros(1, length(new_node_ids));
            count = count + 1;
            loc_set_z_vec = zeros(1, length(nodes_in_range));
            %             loc_set_z_vec = [];
            for k = 1 : length(probable_loc)
                loc = probable_loc(k);
                loc_idx = find(new_node_ids == loc);
                distance = sqrt((new_df_nodes{loc_idx, 2} - new_df_nodes{idx, 2})^2 + (new_df_nodes{loc_idx, 3} - new_df_nodes{idx, 3})^2);
                %                 loc_set_z_vec(1, k) = exp(-distance*100);     %EPSILON = 100
                % loc_set_z_vec(1, loc_idx) = exp(-distance*10);     %EPSILON = 5
                loc_set_z_vec(1, k) = exp(-distance*5);     %EPSILON = 5
            end
            loc_set_z_sum = sum(loc_set_z_vec(1, :));
            loc_set_z_vec(1, :) = loc_set_z_vec(1, :)/loc_set_z_sum;
            z_vec = loc_set_z_vec;
            node_obf_mat = z_vec;
            node_ids_in_range_mat = find(nodes_in_range == id);
        end

        % obfuscation of node
        obf_node = 1;
        obf_node = obfuscation_node(new_df_nodes, z_vec);      
        obf_traj{end+1} = new_node_ids(obf_node);

        % Utility Loss Calculation
        node_utility_loss = utility_loss(node_ids_in_range_mat, :);
        expected_utility_loss = nansum(z_vec .* node_utility_loss);
        expected_utility_loss_traj = [expected_utility_loss_traj, expected_utility_loss];

        % Posterior Computation
        denomi = [];
        % Real_idx = find(ismember(new_node_ids, nodes_in_range));
        obf_loc_posterior = zeros(1, length(new_node_ids));
        for k = 1:length(Real_idx)
            R = Real_idx(k);                    % R = Real Location
            obf_prob = z_vec(k);                        % obf_prob = Probability from obfuscation matrix for OI in RLV
            priorDistribution_R = prior_distribution(R);                    % prior_distribution_R = prior distribution for k_th node in R
            denomi = [denomi, (obf_prob * priorDistribution_R)];                                  % Bayes formula's denominator
        end
        for l = 1:length(Real_idx)
            obf_prob = z_vec(l);
            R = Real_idx(l);                    % R = Real Location
%             r_idx = find(new_node_ids == R);
            priorDistribution_R = prior_distribution(R);
            posterior_prob = (obf_prob * priorDistribution_R)/sum(denomi);                                   % Bayes Formula
            obf_loc_posterior(1, R) = posterior_prob;
        end
        traj_posts{end+1} = obf_loc_posterior;
% 
    end
    obfuscated_traj{end+1,1} = obf_traj;
    obf_mat_traj{end+1} = node_obf_mat;
    utility_loss_mat{end+1} = node_utility_loss_mat;
    expected_utility_loss_alltrajs{end+1, 1} = expected_utility_loss_traj;
    average_expected_utility_loss= [average_expected_utility_loss, nanmean(expected_utility_loss_traj)];
    post_vectors{end+1, 1} = traj_posts;
end
function d = haversine_calc(lat1, lon1, lat2, lon2)
    R = 6372.8; % Earth's radius in kilometers
    dLat = deg2rad(lat2 - lat1);
    dLon = deg2rad(lon2 - lon1);
    lat1 = deg2rad(lat1);
    lat2 = deg2rad(lat2);
    a = sin(dLat/2)^2 + sin(dLon/2)^2*cos(lat1)*cos(lat2);
    c = 2*asin(sqrt(a));
    d = R*c;
end


% 
% adjacencyMatrix = adjacency(G);
% nodes_num = [];
% range_nodes = {};
% zeroNodes = [];
% for i = 1 : length(new_node_ids)
%     
%     adjacencyMatrix = adjacency(G);
%     i = 1;
%     disp(i)
%     centralNode = new_node_ids(i);
%     radius = 3;
% 
%     % Example: Compute distances using breadth-first search
%     distancesToCentralNode = graphshortestpath(adjacencyMatrix, centralNode);
% 
%     % Example: Find nodes within the specified radius
%     nodesWithinRadius = find(distancesToCentralNode <= radius);
%     nodesWithinRadius_idx = find(ismember(new_node_ids, nodesWithinRadius));
%     figure;
%     scatter(new_df_nodes.x, new_df_nodes.y, 'filled');
%     hold on;
%     scatter(new_df_nodes.x(i), new_df_nodes.y(i), 'r', 'filled');
%     hold on;
%     % Highlight nodes within the specified radius
%     scatter(new_df_nodes.x(nodesWithinRadius_idx), new_df_nodes.y(nodesWithinRadius_idx), 'g', 'filled');
% 
%     if(length(nodesWithinRadius) < 2)
%         zeroNodes = [zeroNodes, centralNode];
%     end
%     nodes_num = [nodes_num, length(nodesWithinRadius)];
%     range_nodes{end+1, 1} = nodesWithinRadius;
% 
% end
% 
% % Check if there is any value less than 5
% isValueLessThan5 = any(nodes_num < 2);
% 
% % Display the result
% if isValueLessThan5
%     disp('There is at least one value less than 5 in the array.');
% else
%     disp('All values are greater than or equal to 5 in the array.');
% end





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
%         obf_traj{end+1} = new_node_ids(obf_node);
% 
%     end
%     obfuscated_traj{end+1,1} = obf_traj;
% end




% obf_x = [];
% obf_y = [];
% real_x = [];
% real_y = [];
% for i = 1 : length(new_traj{2,1})
%     real_node_id = new_traj{2,1}(i);
%     idx = find(new_node_ids == real_node_id);
%     real_x = [real_x, new_df_nodes{idx, "x"}];
%     real_y = [real_y, new_df_nodes{idx, "y"}];
% 
%     obf_node_id = obfuscated_traj{1,1}(i);
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