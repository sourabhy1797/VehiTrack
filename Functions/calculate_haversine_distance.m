function [trajectories, traject_coord] = calculate_haversine_distance(dataset_trajectories, col_x, col_y)
%% Objective:
    % This function is to approximate the obtained x and y coordinates from dataset according to the openstreet map and then
    % converts the trajectories in that dataset into particular format
    % which can used for other computations. This method will convert the
    % dataset into two column structure where first column contains the
    % vectors of trajectories, and another column contains the time stamp.
%% Input:
    % dataset_trajectories: These are mainframe trajectories recorded from the main taxi_sample_2.mat dataset
    % col_x: It contains all the x co-ordinates of the all the nodes retrieved from openstreet api
    % col_y: It contains all the y co-ordinates of the all the nodes retrieved from openstreet api
%% Output:
    % trajectories: Matrix which contains trajectory nodes in first column and timestamps in second column
    % traject_coord: Matrix which contains x and y coordinates of trajectory nodes
    
    addpath('./Functions/haversine/'); 
    trajectories = {};
    traject_coord = {};
    for i = 1:length(dataset_trajectories)
        i
        % Initialize list to store closest nodes for each selected node
        traj_x = dataset_trajectories(i).x;
        traj_y = dataset_trajectories(i).y;
        closest_nodes = zeros(size(traj_x));
        closest_nodes_x = zeros(size(traj_x));
        closest_nodes_y = zeros(size(traj_x));

        % Loop through each selected node
        for j = 1:length(traj_x)
            selected_node_x = traj_x(j);
            selected_node_y = traj_y(j);
            min_dist = inf; % Initialize minimum distance to infinity
            closest_node = -1; % Initialize closest node to -1 (invalid node index)
            closest_node_x = -1;
            closest_node_y = -1;

            % Loop through all nodes and calculate Haversine distance to selected node
            for k = 1:length(col_x)
                all_node_x = col_x(k);
                all_node_y = col_y(k);
                distance = haversine_calc(selected_node_x, selected_node_y, all_node_x, all_node_y);
                if distance < min_dist % Check if distance is smaller than current minimum
                    min_dist = distance;
                    closest_node = k; % Update closest node index
                    closest_node_x = col_x(k);
                    closest_node_y = col_y(k);
                end
            end

            closest_nodes(j) = closest_node; % Store closest node index for selected node
            closest_nodes_x(j) = closest_node_x; % Store closest node x for selected node
            closest_nodes_y(j) = closest_node_y; % Store closest node y for selected node
        end
        
        trajectories{i,1} = closest_nodes;
        traject_coord{i,1} = closest_nodes_x;
        traject_coord{i,2} = closest_nodes_y;
        
        % Time in base dataset represents traveltime from one node to another,
        % whereas we are considering to timestamp in continuous form, below
        % code will convert them into continuous form.
        testtime = zeros(1, length(dataset_trajectories(i).time) ,'int64');
        for j = 2:length(dataset_trajectories(i).time)
            dt = dataset_trajectories(i).time(j) - dataset_trajectories(i).time(j-1);
            testtime(j) = testtime(j-1) + dt;
        end
        trajectories{i,2} = testtime;

        
            
    end
    
end

