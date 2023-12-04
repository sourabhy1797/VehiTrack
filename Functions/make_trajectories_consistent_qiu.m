function [new_trajectories, new_traject_coord] = make_trajectories_consistent_qiu(trajectories, traject_coord, time_interval)
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
    time_min = 99999; 
    time_max = 0; 
    for i = 1:1:size(trajectories, 1)
        if min(trajectories{1,2}) < time_min
            time_min = min(trajectories{i,2});
        end
        if max(trajectories{1,2}) > time_max
            time_max = max(trajectories{i,2});
        end
    end

    % Count the number of time slots
    NR_SLOTS = ceil(time_max/time_interval); 
    NR_TRAJ = size(trajectories, 1);
    NR_LOC = 43160;

    distribution = sparse(NR_TRAJ, NR_SLOTS); 
    for i = 1:1:size(trajectories, 1)
        
        for j = 1:1:size(trajectories{i,2}, 2)
           
            time = ceil((trajectories{i,2}(j)+1)/time_interval); 
            loc = trajectories{i,1}(j); 
             % [i time]
            distribution(i, time) = loc;
            
        end
    end



end