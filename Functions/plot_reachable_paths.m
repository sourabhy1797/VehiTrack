function p = plot_reachable_paths(G, paths)
%% Objective:
    % This function is to plot all the reachable paths.
%% Input:
    % G: Graph which is going to base of the plot
    % paths: It is vector of vectors which contains all the reachable paths
%% Output:
    % p = final plot in which all the reachable paths are highlighted

    % Plot the graph
    figure;
    p = plot(G);
    % p.NodeLabel = {}; % Remove node labels for clarity

    % Highlight the paths
    for i = 1:length(paths)
%         highlight(p, paths{i}, 'EdgeColor', 'r', 'LineWidth', 2);
        highlight(p, paths{i}, 'NodeColor', 'r', 'MarkerSize', 2);
    end
    
    % Highlight starting nodes of paths
    start_nodes = cellfun(@(x) x(1), paths);
    highlight(p, start_nodes, 'NodeColor', 'g', 'MarkerSize', 10);
end
