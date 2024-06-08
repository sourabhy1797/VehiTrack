function [reachable_ids, shortDist, shortPaths] = get_reachable_nodes(G, v, t_start, t_end)
%% Objective:
    % This function is to deduce the all reachable nodes in the given time
    % frame
%% Input:
    % G: weighted adjacency matrix of the graph
    % v: starting node
    % t_start: start time
    % t_end: end time
%% Output:
    % reachable_ids: Array indicating the nodes reachable from v within the time frame [t_start, t_end]
    % shortDist: Array of travel time reaching from v to all the reachable nodes
    % shortPaths: Array of shortest path from v to reach all the reachable nodes
    
    n = numnodes(G);
    dist = inf(1,n); % initialize distance vector to infinity
    path = cell(1,n);
    dist(v) = 0; % set distance of starting node to 0
    for i = 1:n
        [shortest_path, shortest_dist] = shortestpath(G, v, i);
        dist(i) = shortest_dist;
        path{i} = shortest_path;
    end

    threshold = t_end - t_start; % example threshold value

    % find indices of distances less than threshold
    indices = find(dist < threshold);
    reachable_ids = indices;
    % initialize cell array to store paths with distances less than threshold
    shortPaths = cell(size(indices));
    shortDist = zeros(size(indices));
    % iterate over indices and store corresponding path in shortPaths
    for i = 1:length(indices)
        idx = indices(i);
        shortPaths{i} = path{idx};
    end

    for i = 1:length(indices)
        idx = indices(i);
        shortDist(i) = dist(idx);
    end

end

