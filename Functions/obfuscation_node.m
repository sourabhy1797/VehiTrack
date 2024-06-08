function obf_location = obfuscation_node(df_nodes, z_vec)
%% Description:
    % The obfuscation method is to obfuscate all the nodes deduced from
    % taxi sample data i.e. trajactory nodes
%% Input
    % df_nodes: a MATLAB table representing a list of nodes with actual x and y 
    % coordinates and unique osmid; the dataset from openstreet map
    % node_for_obf: a Node which needs to be obfuscated

%% Output
    % all_obf_loc: All the obfuscated location list for all the deduced trajectory ids
    % all_obf_vec: Obfuscation Matrix
    % reachable_vec: All probable reachable nodes within range in whole big ROME Network
%%    
    EPSILON = 100;
    % load('obf_mat_e50.mat')

    % Coordinate matrix from openstreet dataset
    id_x = df_nodes.x;
    id_y = df_nodes.y;
    coordMatrix = horzcat(id_x, id_y);
    NR_LOC = length(id_x);
%     [z_vec, reachable_vec] = obfmatrix_generator_Laplace(coordMatrix, node_for_obf, EPSILON, NR_LOC);
%     z_vec = z_vector(node_for_obf, :);
%     all_obf_loc = [];
%     for t = 1:1:20
        obf_vec = z_vec;
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
                obf_location = (k);
            end
        end
%         all_obf_loc = horzcat(all_obf_loc, obf_location);
%     end
    
end