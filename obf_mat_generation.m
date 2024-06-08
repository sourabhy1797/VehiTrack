%% This file is for geneerating the obfuscation matrix

addpath('./Functions/');

opts = detectImportOptions('./Dataset/smaller_nodes_rome.csv');
opts = setvartype(opts, 'osmid', 'int64');
opts = setvartype(opts, 'x', 'double');
opts = setvartype(opts, 'y', 'double');
df_nodes = readtable('./Dataset/smaller_nodes_rome.csv', opts);
df_nodes_big = readtable('./Dataset/nodes_rome.csv', opts);
df_edges = readtable('./Dataset/edges_rome.csv');

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

id_x = matchingNodes.x;
id_y = matchingNodes.y;

% id_x = new_df.x;
% id_y = new_df.y;

coordMatrix = horzcat(id_x, id_y);




coordinate = coordMatrix;
EPSILON = 7.5;
NR_LOC = length(id_x);
OBF_RANGE = 0.01;
for i = 1:1:NR_LOC
    disp(i)
    for j = 1:1:NR_LOC
        distance = sqrt((coordinate(i, 1) - coordinate(j, 1))^2 + (coordinate(i, 2) - coordinate(j, 2))^2);
        % distance = haversine_calc(coordinate(i, 1), coordinate(i, 2), coordinate(j, 1), coordinate(j, 2));
        if distance <= OBF_RANGE
            z_vector(i, j) = exp(-distance*EPSILON);        % changed i to 1
            reachable_vec(i,j) = 1;
        else
            z_vector(i, j) = 0;             % changed i to 1
            reachable_vec(i,j) = 0;
        end
    end
    z_sum = sum(z_vector(i, :));            %changed i to 1
    z_vector(i, :) = z_vector(i, :)/z_sum;  %changed i to 1
    


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