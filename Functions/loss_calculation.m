function loss = loss_calculation(df_nodes, new_id_list, all_obf_loc, col_x, col_y)
%% Description:
    % The loss_calculation function is to calculate the loss between obfuscated
    % nodes and ccorresponding real nodes
%% Input
    % df_nodes: a MATLAB table representing a list of nodes with actual x and y 
    % coordinates and unique osmid; the dataset from openstreet map
    % new_id_list: a MATLAB array of unique osmid values corresponding to x and y coordinates with duplicates removed
    % all_obf_loc: Corresponding obfuscated location for ids
    % col_x: all x coordinates of nodes of graph G
    % col_y: all y coordinates of nodes of graph G

%% Output
    % loss: Array which contains the loss between actual node and its corresponding obfuscated node
%%
    id_list = df_nodes.osmid;
    id_list_index = zeros(size(id_list));
    for i = 1:1:size(id_list, 1)
        id_list_index(i, 1) = i;
    end
    [~, new_trajectory_ids] = ismember(new_id_list, id_list_index);
    new_trajectory_ids_t = new_trajectory_ids';
    sample_location = 1;
    x_sample = col_x(sample_location);
    y_sample = col_y(sample_location);
    % nodes_distances = distances(G, id_list_index, id_list_index);
    % nodes_distances = distances(G, sample_location, obf_vec(j,1));
    % loss = zeros(length(new_trajectory_ids), length(all_obf_loc(1, :)));
    loss = [];
    for i = 1:1:length(all_obf_loc(1, :))
        loss_obf = zeros(1, length(new_trajectory_ids));
        obf_vec = all_obf_loc(:, i);
        act_vec = new_trajectory_ids_t;
        x_obf = col_x(obf_vec);
        x_act = col_x(act_vec);
        y_obf = col_y(obf_vec);
        y_act = col_y(act_vec);
        for j = 1:length(obf_vec)
    %         loss_obf(1,j) = abs(nodes_distances(sample_location, obf_vec(j,1)) - nodes_distances(sample_location, act_vec(j,1)));
            R = 6371;
            % Define the latitude and longitude coordinates of the two nodes
            lat1 = x_obf(j); % node 1 latitude
            lon1 = y_obf(j); % node 1 longitude
            lat2 = x_act(j); % node 2 latitude
            lon2 = y_act(j); % node 2 longitude
            % Calculate the haversine distance between the two points
            d1 = 2 * R * asin(sqrt(sin(deg2rad(lat2-x_sample)/2)^2 + cos(deg2rad(x_sample)) * cos(deg2rad(lat2)) * sin(deg2rad(lon2-y_sample)/2)^2));
            d2 = 2 * R * asin(sqrt(sin(deg2rad(lat1-x_sample)/2)^2 + cos(deg2rad(x_sample)) * cos(deg2rad(lat1)) * sin(deg2rad(lon1-y_sample)/2)^2));
            loss_obf(1,j) = abs(d1 - d2);
%             loss_obf(1,j) = abs(distances(G, sample_location, obf_vec(j,1)) - distances(G, sample_location, act_vec(j,1)));
        end
        loss = horzcat(loss, loss_obf');
    end
end