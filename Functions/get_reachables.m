function [reachable_ids, reachable_ids_x, reachable_ids_y, shortDist] = get_reachables(G, v, v_x, v_y, t_start, t_end, col_x, col_y)
%% Objective:
    % This function is to deduce the all reachable nodes in the given time
    % frame
%% Input:
    % G: weighted adjacency matrix of the graph
    % v: starting node
    % v_x: x coordinate of the v node
    % v_y: y coordinate of v node
    % t_start: start time
    % t_end: end time
    % col_x: all x coordinates of nodes of graph G
    % col_y: all y coordinates of nodes of graph G
%% Output:
    % reachable_ids: Array indicating the nodes reachable from v within the time frame [t_start, t_end]
    % reachable_ids_x: Array indicating the x coordinates of nodes reachable from v within the time frame [t_start, t_end]
    % reachable_ids_y: Array indicating the y coordinates of nodes reachable from v within the time frame [t_start, t_end]
    % shortDist: Array of travel time reaching from v to all the reachable nodes
    
    distance_threshold = 10; 

    n = numnodes(G);
    nodes = [];
    nodes_x = [];
    nodes_y = [];
    dist = []; % initialize distance vector to infinity
    path = {};
    % assume that G is the graph object with x and y coordinates stored as node attributes
    radius = 0.5; % radius of the circle in km


    center = [v_x, v_y];
    center_rad = deg2rad(center);

    R = 6371; % km
    x_rad = deg2rad(col_x);
    y_rad = deg2rad(col_y);
    delta_x = x_rad - deg2rad(center(1));
    delta_y = y_rad - deg2rad(center(2));
    a = sin(delta_y/2).^2 + cos(center_rad(2)) .* cos(y_rad) .* sin(delta_x/2).^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    d = R * c;
    idx = find(d <= radius);

    
    % find all the nodes within the circle using rangesearch
%     [idx,~] = rangesearch([col_x col_y],[x0 y0],radius_deg, 'Distance',@haversine_calc);
    
    % idx will contain a cell array where each cell contains the node IDs within the circle
    % you can convert this to a single vector using cell2mat
    nodesInCircle = (idx);
    nodesInCircle_x = col_x(nodesInCircle);
    nodesInCircle_y = col_y(nodesInCircle);
    
%     disp(length(nodesInCircle));
    
    
    
    for i = 1:length(nodesInCircle)
        [shortest_path, shortest_dist] = shortestpath(G, v, nodesInCircle(i));
        if shortest_dist < distance_threshold
            nodes = [nodes nodesInCircle(i)];
            nodes_x = [nodes_x nodesInCircle_x(i)];
            nodes_y = [nodes_y nodesInCircle_y(i)];
            dist = [dist shortest_dist];
            path = {path shortest_path};
        end

    end

    threshold = t_end - t_start; % example threshold value
    % find indices of distances less than threshold
    indices = find(dist < threshold);
    reachable_ids = [];
    reachable_ids_x = [];
    reachable_ids_y = [];
    shortDist = [];
    for i = 1:length(indices)
        reachable_ids = [reachable_ids nodes(indices(i))];
        reachable_ids_x = [reachable_ids_x nodes_x(indices(i))];
        reachable_ids_y = [reachable_ids_y nodes_y(indices(i))];
        shortDist = [shortDist dist(indices(i))];
    end

end