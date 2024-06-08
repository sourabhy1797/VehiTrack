function violation_percentages = compute_violation_percentages_LP(obfuscation_matrices, trajectory, target_region_nodes, epsilon, distances)
    num_nodes = length(target_region_nodes); % Number of nodes in the target region
    num_trajectory_nodes = length(trajectory);
    
    violation_percentages = zeros(num_trajectory_nodes, 1);
    
    for i = 1:num_trajectory_nodes
        current_node = trajectory(i);
        current_target_region_nodes = target_region_nodes{i};
        num_current_target_nodes = length(current_target_region_nodes);
        num_non_target_nodes = num_nodes - num_current_target_nodes;
        
        % Access obfuscation matrix for the current node
        obfuscation_matrix = obfuscation_matrices{i};
        
        non_target_probability_sum = sum(obfuscation_matrix(current_node, :));
        target_probability_sum = sum(obfuscation_matrix(:, current_node));
        
        % Calculate violation for each target region node
        for j = 1:num_current_target_nodes
            target_node = current_target_region_nodes(j);
            distance = distances(current_node, target_node);  % Distance between current node and target node
            
            if target_probability_sum > 0
                expected_probability = target_probability_sum * exp(distance * epsilon);
                if non_target_probability_sum > expected_probability
                    violation_percentages(i) = violation_percentages(i) + (non_target_probability_sum - expected_probability) / num_non_target_nodes;
                end
            end
        end
        
        violation_percentages(i) = violation_percentages(i) * 100;  % Convert to percentage
    end
end