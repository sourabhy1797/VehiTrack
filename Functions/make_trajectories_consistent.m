function [new_trajectories, traject_coord] = make_trajectories_consistent(trajectories, traject_coord, G, col_x, col_y, time_interval)
%% Description:
    % The preprocessing_file function is defined to convert the
    % base trajectories into consistent trajectories based on openstreet
    % map based graph considering time stamps as major field or
    % charaterstics.
%% Input
    % trajectories: Old trajectories IDs and time stamps which needs to be
                % preprocessed
    % traject_coord: Old trajectories x and y coordinates which needed to
                % be preprocessed
    % G: Graph computed using openstreet map
%% Output
    % new_trajectories: These are the trajectories which are preprocessed
                    % form of base datasset trajectories, and these tajectories are
                    % consistent in terms of time.
    % traject_coord: Matrix which contains x and y coordinates of trajectory nodes
    
    %% [?? Chenxi] Make sure it covers the larger interval that needs to be filled by multiple locations. 
    % create a copy of the trajectories
    new_trajectories = trajectories;
    disp(new_trajectories)

    % loop through each trajectory
    for j = 1:length(new_trajectories)
        j

        % extract the trajectory, time stamps and x, y coordinates from the current trajectory
        t = double(new_trajectories{j,2});
        traj = new_trajectories{j,1};
        traj_x = traject_coord{j,1};
        traj_y = traject_coord{j,2};

        % initialize new time stamps and empty new trajectory
        t_new = [0];
        new_z = [];

        % loop through each time stamp and interpolate between them
        for i = 2: length(t)

            % create an array with the start and end time stamp
            arr = [t(i-1) t(i)];

            % calculate the time difference between start and end
            arr1 = [t(i-1)-t(i-1) t(i)-t(i-1)];

            % calculate how many 10-second intervals are between the start and end time stamp
            num = (floor(arr1(2)/time_interval));
            if num>0

                % if there are multiple 10-second intervals, insert new time stamps
                ins = [];
                for k = 1:num
                    ins = [ins (k*time_interval)];
                end
                arr2 = [ins arr1(end)];
                arr2 = arr2 + t(i-1);
                t_new = [t_new arr2];
            else
                % if there is only one time stamp, append it to the new time stamps
                t_new = [t_new arr(2)];
            end

            % remove any duplicates from the new time stamps
            t_new = (unique(t_new));
        end

        % map the new time stamps to their corresponding index in the original time stamp vector
        [~, idx] = ismember(t_new, t);
        
        % create a new trajectory and x, y coordinates based on the new time stamps
        z_new = idx;
        traj_x_new = idx;
        traj_y_new = idx;
        insert_idx = [];

        % loop through each index in the new time stamp vector and populate the new trajectory and x, y coordinates
        for i = 1 : length(z_new)
            if(z_new(i)>0)
                ind = z_new(i);
                z_new(i) = traj(ind);
                traj_x_new(i) = traj_x(ind);
                traj_y_new(i) = traj_y(ind);
            elseif z_new(i) == 0
                insert_idx = [insert_idx i];
            end
        end
        traj = z_new;
        t = t_new;
        traj_x = traj_x_new;
        traj_y = traj_y_new;

        % interpolate between missing time stamps
        for i = 1:length(insert_idx)
            t_end = t(insert_idx(i));
            t_start = t(insert_idx(i)-1);
            n = traj(insert_idx(i)-1);
            n_x = traj_x(insert_idx(i)-1);
            n_y = traj_y(insert_idx(i)-1);
            % deducing the reachable nodes within that time frame
            [reachable_nodes, reachable_nodes_x, reachable_nodes_y, shortDist] = get_reachables(G, n, n_x, n_y, t_start, t_end, col_x, col_y);    
            shortDistNext = [];
            idx = find(shortDist < (t_end - t_start), 1, 'last');

            % insert the new trajectory and x, y coordinates
            traj(insert_idx(i)) = reachable_nodes(idx);
            traj_x(insert_idx(i)) = reachable_nodes_x(idx);
            traj_y(insert_idx(i)) = reachable_nodes_y(idx);

        end

        % Get the maximum value of the time stamps and calculate the number of intervals (k_lim) with a length of 10.
        k_lim = floor(t(end)/time_interval);

        % Initialize new_t and t_idx arrays.
        new_t = [];
        t_idx = [];
        k = 0;

        % Initialize the median array to 5.
        med_arr1 = time_interval/2;

        % Loop through each interval.
        while k<=k_lim

            % Create an array (arr) containing all time stamps within the current interval.
            arr = [];
            arr_idx = [];
            for i = 1:length(t)
                if(floor(t(i)/time_interval) == k)
                    arr = [arr t(i)];
                    arr_idx = [arr_idx i];
                end
            end

            % Calculate the value of the median and its corresponding index in the original array.
            diff = abs(arr - med_arr1);
            [~, idx] = min(diff);
            med_arr = arr(idx);
            med_idx = find(arr == med_arr);
            med_idx = arr_idx(med_idx);

            % Add the median value and index to the new_t and t_idx arrays.
            k = k+1;
            med_arr1 = med_arr1 + time_interval;
            new_t = [new_t med_arr];
            t_idx = [t_idx med_idx];
        end

        % Update the trajectory array with the values at the selected indices.
        updated_traj = traj(t_idx);
        updated_traj_x = traj_x(t_idx);
        updated_traj_y = traj_y(t_idx);

        %% [?? Chenxi] What if an interval includes multiple locations? Does this code cover the case? 
        new_trajectories{j,1} = updated_traj;
        new_trajectories{j,2} = new_t;
        traject_coord{j,1} = updated_traj_x;
        traject_coord{j,2} = updated_traj_y;
    end
end