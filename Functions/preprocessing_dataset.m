function [trajectories, traject_coord] = preprocessing_dataset(df_nodes, vehicle_trace)
%% Description:
    % The preprocessing_dataset function read the given taxi dataset, 
    % approximate the obtained x and y coordinates from dataset according to the openstreet map and then
    % converts the trajectories in that dataset into particular format
    % which can used for other computations. This method will convert the
    % dataset into two column structure where first column contains the
    % vectors of trajectories, and another column contains the time stamp.
%% Input
    % df_nodes: Openstreet map nodes details
    % vehicle_trace: Structure containing the trajectory data
    
%% Output
    % trajectories: Matrix which contains trajectory nodes in first column and timestamps in second column
    % traject_coord: Matrix which contains x and y coordinates of trajectory nodes
    
%% Reading and preprocessing the dataset file and dataset in it
%     load('.\Dataset\trace_taxi_sample_2.mat');
    trajectory_table= {};
    for i = 1:length(vehicle_trace)
        i
        df_trace_taxi_sample = struct2table(vehicle_trace(i));
        df_trace_taxi_sample1 = table2array(df_trace_taxi_sample);
        column1 = df_trace_taxi_sample1(:,1);
        column3 = df_trace_taxi_sample1(:,3);
        column4 = df_trace_taxi_sample1(:,4);
        trajectory_table{i,1} = i;
        trajectory_table{i,2} = column4';
        trajectory_table{i,3} = column3';
        trajectory_table{i,4} = column1';
    end
    myFieldNames = {'ID', 'x', 'y', 'time'};

    % create a table with the data and field names
    myTable = cell2table(trajectory_table, 'VariableNames', myFieldNames);

    % convert the table to a structure
    dataset_trajectories = table2struct(myTable);


    %% Section to deduce the node ids according openstreet map by performing Approximation


    % Extract the column as an array
    col_x = table2array(df_nodes(:, 'x'));  % Actual x coordinate from the nodes data
    col_y = table2array(df_nodes(:, 'y'));  % Actual y coordinate from the nodes data
    col_osmid = table2array(df_nodes(:, 'osmid'));  % Actual unique osmid from the nodes data

    [trajectories, traject_coord] = calculate_haversine_distance(dataset_trajectories, col_x, col_y);
    
    
    
%% Ignore this portion of code. It is the old approach to calculate the approximate nodes. It uses euclidean distance for calculating the approximate node.
%     trajectories = {};
%     traject_coord = {};
    
    

%     for i = 1:length(dataset_trajectories)
%         closest_ids = zeros(1, length(dataset_trajectories(i).x) ,'int64');
%         closest_x = zeros(1, length(dataset_trajectories(i).x));
%         closest_y = zeros(1, length(dataset_trajectories(i).x));
%         for j = 1:length(dataset_trajectories(i).x)
%             % Calculate distances between actual x and y co-ordinates and all assumed x and y co-ordinates
%             dx = dataset_trajectories(i).x(j) - col_x;
%             dy = dataset_trajectories(i).y(j) - col_y;
%             distances = sqrt(dx.^2 + dy.^2);
% 
%             % Find index of closest assumed x and y co-ordinates
%             [~, idx] = min(distances);
%             closest_ids(j) = int64((idx));
%             closest_x(j) = (col_x(idx));
%             closest_y(j) = (col_y(idx));
%         end
%         trajectories{i,1} = closest_ids;
%         traject_coord{i,1} = closest_x;
%         traject_coord{i,2} = closest_y;
% 
%         % Time in base dataset represents traveltime from one node to another,
%         % whereas we are considering to timestamp in continuous form, below
%         % code will convert them into continuous form.
%         testtime = zeros(1, length(dataset_trajectories(i).time) ,'int64');
%         for j = 2:length(dataset_trajectories(i).time)
%             dt = dataset_trajectories(i).time(j) - dataset_trajectories(i).time(j-1);
%             testtime(j) = testtime(j-1) + dt;
%         end
%         trajectories{i,2} = testtime;
%     end
end