function[new_x_list, new_y_list, new_id_list] = approximation(df, df_nodes)
%% Description:
    % The approximation function takes in a matrix of data and a list of nodes 
    % with actual x and y coordinates, and calculates the closest assumed x and
    % y coordinates for each point in the data matrix. It then removes duplicates 
    % and returns the final lists of x and y coordinates with corresponding unique osmid values.
%% Input
    % df: a MATLAB table representing a matrix of data which is dataset we
    % been given of taxi
    % df_nodes: a MATLAB table representing a list of nodes with actual x and y 
    % coordinates and unique osmid; the dataset from openstreet map

%% Output
    % new_x_list: a MATLAB array of x coordinates with duplicates removed
    % new_y_list: a MATLAB array of y coordinates with duplicates removed
    % new_id_list: a MATLAB array of unique osmid values corresponding to x and y coordinates with duplicates removed
%%    
    % Convert table to matrix
    M = table2array(df);
    M(isnan(M)) = 999;
    
    % Extract the column as an array
    col_x = table2array(df_nodes(:, 'x'));  % Actual x coordinate from the nodes data
    col_y = table2array(df_nodes(:, 'y'));  % Actual y coordinate from the nodes data
    col_osmid = table2array(df_nodes(:, 'osmid'));  % Actual unique osmid from the nodes data
    
    % Convert the array to a list
    actual_x = cell2mat(num2cell(col_x));
    actual_y = cell2mat(num2cell(col_y));
    osmid = cell2mat(num2cell(col_osmid));
    
    % Creating seperate empty list for storing edges x and y co-ordinates
    % for start point u and end point v
    ux = [];
    vx = [];
    uy = [];
    vy = [];
    
    for i = 1:size(M,1)
        testux = [];
        testvx = [];
        testuy = [];
        testvy = [];
        for j = 2:64            % Till 64 all are x values and if the value is not null then that value is stored in normal testux list
            if M(i,j) ~= 999
               testvx = horzcat(testvx, M(i,j));
               testux = horzcat(testux, M(i,j));
            end
        end
        % Each row in that dataset have both x 
        % and y values, let take x values, suppose there are 3 values and rest 61 is null.
        % There 1 will have edge with 2 and 2 will have edge 3 but 3 cannot be initiating point for any edge for that data set.
        % Therefore ux till second last element, and vx from second to last element
        ux = [ux, testux(1:end-1)]; 
        vx = [vx, testvx(2:end)];
    
        % Same concept for y
        for j = 65:127
            if M(i,j) ~= 999
               testvy = horzcat(testvy, M(i,j));
               testuy = horzcat(testuy, M(i,j));
            end
        end
        uy = [uy, testuy(1:end-1)];
        vy = [vy, testvy(2:end)];
    end
    
    % Pre-allocate closest assumed x and y co-ordinates
    closest_x = zeros(length(ux), 1);
    closest_y = zeros(length(uy), 1);
    closest_osmid = zeros(length(ux), 1, 'int64');
    % Loop over actual x and y co-ordinates
    for i = 1:length(ux)
        % Calculate distances between actual x and y co-ordinates and all assumed x and y co-ordinates
        dx = ux(i) - actual_x;
        dy = uy(i) - actual_y;
        distances = sqrt(dx.^2 + dy.^2);
        
        % Find index of closest assumed x and y co-ordinates
        [~, idx] = min(distances);
        
        % Assign closest assumed x and y co-ordinates to corresponding index
        closest_x(i) = actual_x(idx);
        closest_y(i) = actual_y(idx);
        closest_osmid(i) = (osmid(idx));
    end
    
    
    
    % To drop the repeatative values
    % Initialize hash table
    hash_table = containers.Map();
    
    % Initialize new lists
    new_x_list = [];
    new_y_list = [];
    new_id_list = int64([]);
    
    % Loop through lists and add non-duplicates to new lists
    for i = 1:length(closest_x)
        x = closest_x(i);
        y = closest_y(i);
        id = closest_osmid(i);
        key = sprintf('%d_%d', x, y); % create hash table key
        
        if ~isKey(hash_table, key)
            % add to hash table and new lists
            hash_table(key) = true;
            new_x_list(end+1) = x;
            new_y_list(end+1) = y;
            new_id_list(end+1) = id;
        end
    end
end
