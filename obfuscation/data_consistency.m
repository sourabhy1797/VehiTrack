%% This file is to bring data consistency

addpath('./Functions/');

opts = detectImportOptions('./Dataset/smaller_nodes_rome.csv');
opts = setvartype(opts, 'osmid', 'int64');
opts = setvartype(opts, 'x', 'double');
opts = setvartype(opts, 'y', 'double');
df_nodes = readtable('./Dataset/smaller_nodes_rome.csv', opts);
% df_nodes_big = readtable('./Dataset/nodes_rome.csv', opts);
df_edges = readtable('./Dataset/edges_rome.csv');

%% New Nodes Dataset
new_df = readtable('./Dataset/node.csv');
newColumnOrder = {'node', 'lat', 'lng'};
new_df = new_df(:,newColumnOrder);
new_df = renamevars(new_df, 'lat', 'y');
new_df = renamevars(new_df, 'lng', 'x');

columnsToKeep = {'y', 'x'};

% Extract the relevant columns (Node_ID, Latitude, Longitude)
% nodes_data_bigger = df_nodes_big(:, columnsToKeep);

% id_list = [];
% for i = 1:height(nodes_data_bigger)
%     id_list = [id_list, i];
% end
% nodes_data_bigger = addvars(nodes_data_bigger, id_list', 'NewVariableName', 'Node_ID');

nodes_data_smaller = df_nodes(:, columnsToKeep);

% large_coord = df_nodes_big(:, 2:3);
large_coord = new_df(:, 2:3);           % here used the new nodes dataset
smaller_coord = df_nodes(:,2:3);

% % Use ismember to find the indices of matching coordinates
[~, indices] = ismember(large_coord, smaller_coord, 'rows');

% Create a new dataset containing the matching nodes
matchingNodes = new_df(indices > 0, :);
tab2arr = table2array(matchingNodes);


%% New Trajectory Dataset
trajectoryDataset = readtable('./Dataset/matching_result.csv');
trajectoryStrings = trajectoryDataset.Node_list;
newTrajectories = cell(length(trajectoryStrings), 1);

for i = 1 : numel(trajectoryStrings)
    trajectoryArray = str2num(strrep(trajectoryStrings{i}, ',', ' '));
    newTrajectories{i,1} = trajectoryArray;
end