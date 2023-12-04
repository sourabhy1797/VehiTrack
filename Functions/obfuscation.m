function [all_obf_loc, all_obf_vec ]= obfuscation(df_nodes, new_id_list)
%% Description:
    % The obfuscation method is to obfuscate all the nodes deduced from
    % taxi sample data i.e. trajactory nodes
%% Input
    % df_nodes: a MATLAB table representing a list of nodes with actual x and y 
    % coordinates and unique osmid; the dataset from openstreet map
    % new_id_list: list of updated trajectory node ids which changed due
    % that indexing thing.

%% Output
    % all_obf_loc: All the obfuscated location list for all the deduced trajectory ids
%%    
    EPSILON = 100;
    id_list = df_nodes.osmid;
    id_list_index = zeros(size(id_list));
    for i = 1:1:size(id_list, 1)
        id_list_index(i, 1) = i;
    end
    new_id_list = new_id_list';
    [~, new_trajectory_ids] = ismember(new_id_list, id_list_index);
    % Since the graph is made using the indices of the node ids, therefore
    % recoridng the node id indices for calculating the node distances.
    
    
    % Coordinate matrix from openstreet dataset
    id_x = df_nodes.x;
    id_y = df_nodes.y;
    coordMatrix = horzcat(id_x, id_y);
    NR_LOC = length(id_x);
    all_obf_vec = [];
    for i = 1:1:length(new_trajectory_ids)
        z_vec = obfmatrix_generator_Laplace(coordMatrix, new_trajectory_ids(i), EPSILON, NR_LOC);
        all_obf_vec = [all_obf_vec; z_vec];
    end
    
    all_obf_loc = [];
    for t = 1:1:10        % Repeative loop for 10 or 100 or n number times the obfuscated location needed for a particular trajectory location
        obf_location = zeros(length(new_trajectory_ids),1);
        for i = 1:1:length(new_trajectory_ids)
            obf_vec = all_obf_vec(i, :);
            % Creating a sum array
            sum_array = zeros(length(obf_vec)+1, 1);
            rs = 0;
            for j = 1:length(obf_vec)
                rs = rs + obf_vec(j);
                sum_array(j+1) = rs;
            end
            % Applying randomization operation to select the interval
            rand_interval = rand(1,1);
            for k = 1:1:(length(sum_array)-1)
                if(rand_interval >= sum_array(k) && rand_interval < sum_array(k+1))
                    obf_location(i) = (k);
                end
            end
        end
        all_obf_loc = horzcat(all_obf_loc, obf_location);
    end
end