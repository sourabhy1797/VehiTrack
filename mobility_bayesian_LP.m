%% This is mainstream file for generating the 1 posteriors for 1 node in trajectory for benchmark

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
df_edges_big = readtable('./Dataset/edges_rome.csv');
df_edges = readtable('./Dataset/smaller_edges_rome.csv');

load('.\Dataset\rome_taxicab\traject_coord_approx.mat');
load('.\Dataset\rome_taxicab\trajectories_approx.mat'); 

load('.\obf_traj_data\obf_traj_e5.mat');


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


nonMatchingNodes = new_nodes_df(indices == 0, :);
nonMatchingNodesIds = (nonMatchingNodes.node);

edgesToRemove = new_edges_df(ismember(new_edges_df.s_node, nonMatchingNodesIds) | ismember(new_edges_df.e_node, nonMatchingNodesIds), :);

new_edges_data = new_edges_df(~ismember(new_edges_df.s_node, nonMatchingNodesIds) & ~ismember(new_edges_df.e_node, nonMatchingNodesIds), :);


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

for i = 1:length(new_traj)
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




%% Calculation of Posteriors for each nodes based its occurence in trajectory
pois = readtable('./Dataset/pois.csv');
pois = table2array(pois(:,1));
new_pois = intersect(new_node_ids, pois);
new_pois_idx = find(ismember(new_node_ids, new_pois));

post_vectors = {};
for i = 1 : 100
    disp(i)
    node_utility_loss_mat = {};
    traj = new_traj{i,1};
    obfuscated_traj = obf_traj{i,1};
    traj_posts = {};
    id = traj(1);
    obf_id = obfuscated_traj{1};
    idx = find(new_node_ids == id);
    probable_loc = [];

    dist = [];
    for k = 1 : length(new_node_ids)
        [path, distances] = shortestpath(G, obf_id, new_node_ids(k));
        dist = [dist, distances];
        if distances < 250
            probable_loc = [probable_loc, new_node_ids(k)];
        end
    end
    elementsToRemove = probable_loc(~ismember(probable_loc, new_node_ids));
    probable_loc = setdiff(probable_loc, elementsToRemove);
    probable_loc = [probable_loc, obf_id];
    nodes_in_range = probable_loc;


    utility_loss = zeros(length(nodes_in_range), length(nodes_in_range));

    % utility loss computation
    for m = 1 : length(nodes_in_range)
        real_id = nodes_in_range(m);
        for n = 1 : length(nodes_in_range)
            obf_id = nodes_in_range(n);
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

    obf_mat = LPObfuscationMatrx(utility_loss, distanceMatrix, 5);
    node_obf_mat{end+1} = obf_mat;
    node_ids_in_range_mat = find(nodes_in_range == id);
    z_vec = obf_mat(node_ids_in_range_mat, :);

    % obfuscation of node
    obf_node = obfuscation_node(new_df_nodes, z_vec);


%     loc_set_z_vec = [];
%     for k = 1 : length(probable_loc)
%         loc = probable_loc(k);
%         loc_idx = find(new_node_ids == loc);
%         distance = sqrt((new_df_nodes{loc_idx, 2} - new_df_nodes{idx, 2})^2 + (new_df_nodes{loc_idx, 3} - new_df_nodes{idx, 3})^2);
%         loc_set_z_vec(1, k) = exp(-distance*10);     %EPSILON = 10
%     end
%     loc_set_z_sum = sum(loc_set_z_vec(1, :));
%     loc_set_z_vec(1, :) = loc_set_z_vec(1, :)/loc_set_z_sum;
    obf_loc_posterior = zeros(1, length(new_node_ids));
    denomi = [];
    for k = 1:length(probable_loc)
        R = probable_loc(k);                    % R = Real Location
        r_idx = find(new_node_ids == R);
        obf_prob = z_vec(k);                        % obf_prob = Probability from obfuscation matrix for OI in RLV
        priorDistribution_R = prior_distribution(r_idx);                    % prior_distribution_R = prior distribution for k_th node in R
        denomi = [denomi, (obf_prob * priorDistribution_R)];                                  % Bayes formula's denominator
    end
    for l = 1:length(probable_loc)
        obf_prob = z_vec(l);
        R = probable_loc(l);                    % R = Real Location
        r_idx = find(new_node_ids == R);
        priorDistribution_R = prior_distribution(r_idx);
        posterior_prob = (obf_prob * priorDistribution_R)/sum(denomi);                                   % Bayes Formula
        obf_loc_posterior(1, r_idx) = posterior_prob;
    end
    % scatter(nodes_x, nodes_y, 30, obf_loc_posterior, "filled");
    traj_posts{end+1} = obf_loc_posterior;

% 
    for j = 2 : length(traj)
        id = traj(j);
        obf_id = obfuscated_traj{j};
        idx = find(new_node_ids == id);
        new_probable_loc = [];
        for k =  1:length(probable_loc)
            for m = 1 : length(new_node_ids)
                [path, distances] = shortestpath(G, probable_loc(k), new_node_ids(m));
                if isfinite(distances)
                    if ~ismember(probable_loc(k), new_probable_loc)
                        [obf_path, obf_dist] = shortestpath(G, probable_loc(k), obf_id);
                        if(obf_dist < 250)
                            new_probable_loc= [new_probable_loc, probable_loc(k)];
                        end
                    end
                end
            end
        end



%         reach_vec = reachable_vec(:, obf_node);
%         Real_idx = find(reach_vec == 1);
% 
%         probable_loc = intersect(Real_idx, new_probable_loc);
        probable_loc = [new_probable_loc,  obf_id];
        nodes_in_range = probable_loc;
        utility_loss = zeros(length(nodes_in_range), length(nodes_in_range));

        % utility loss computation
        for m = 1 : length(nodes_in_range)
            real_id = nodes_in_range(m);
            for n = 1 : length(nodes_in_range)
                obf_id = nodes_in_range(n);
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

        obf_mat = LPObfuscationMatrx(utility_loss, distanceMatrix, 5);
        node_obf_mat{end+1} = obf_mat;
        node_ids_in_range_mat = find(nodes_in_range == id);
        z_vec = obf_mat(node_ids_in_range_mat, :);

        % obfuscation of node
        obf_node = obfuscation_node(new_df_nodes, z_vec);


%         loc_set_z_vec = [];
%         for k = 1 : length(probable_loc)
%             loc = probable_loc(k);
%             loc_idx = find(new_node_ids == loc);
%             distance = sqrt((new_df_nodes{loc_idx, 2} - new_df_nodes{idx, 2})^2 + (new_df_nodes{loc_idx, 3} - new_df_nodes{idx, 3})^2);
%             loc_set_z_vec(1, k) = exp(-distance*10);     %EPSILON = 10
%         end
%         loc_set_z_sum = sum(loc_set_z_vec(1, :));
%         loc_set_z_vec(1, :) = loc_set_z_vec(1, :)/loc_set_z_sum;
        obf_loc_posterior = zeros(1, length(new_node_ids));
        denomi = [];
        for k = 1:length(probable_loc)
            R = probable_loc(k);                    % R = Real Location
            r_idx = find(new_node_ids == R);
            obf_prob = loc_set_z_vec(k);                        % obf_prob = Probability from obfuscation matrix for OI in RLV
            priorDistribution_R = prior_distribution(r_idx);                    % prior_distribution_R = prior distribution for k_th node in R
            denomi = [denomi, (obf_prob * priorDistribution_R)];                                  % Bayes formula's denominator
        end
        for l = 1:length(probable_loc)
            obf_prob = loc_set_z_vec(l);
            R = probable_loc(l);                    % R = Real Location
            r_idx = find(new_node_ids == R);
            priorDistribution_R = prior_distribution(r_idx);
            posterior_prob = (obf_prob * priorDistribution_R)/sum(denomi);                                   % Bayes Formula
            obf_loc_posterior(1, r_idx) = posterior_prob;
        end
        % scatter(nodes_x, nodes_y, 30, obf_loc_posterior, "filled");
        traj_posts{end+1} = obf_loc_posterior;
    end
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
