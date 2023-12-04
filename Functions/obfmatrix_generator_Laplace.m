function [z_vector, reachable_vec] = obfmatrix_generator_Laplace(coordinate, i, EPSILON, NR_LOC)
%% Description:
    % The obfmatrix_generator_laplace is function which generate the
    % obfuscaed location vector
%% Input
    % coordinate: Actual x and y coordinate from the openstreet map data
    % i: location on which obfuscation will happen
    % EPSILON: randomized value
    % NR_LOC: total number of nodes

%% Output
    % z_vector: obfuscation vector for that node
%%
    OBF_RANGE = 100; 
    for j = 1:1:NR_LOC
        distance = sqrt((coordinate(i, 1) - coordinate(j, 1))^2 + (coordinate(i, 2) - coordinate(j, 2))^2); 
        if distance <= OBF_RANGE
            z_vector(1, j) = exp(-distance*EPSILON);        % changed i to 1
            reachable_vec(1,j) = 1;
        else 
            z_vector(1, j) = 0;             % changed i to 1
            reachable_vec(1,j) = 0;
        end
    end
    z_sum = sum(z_vector(1, :));            %changed i to 1
    z_vector(1, :) = z_vector(1, :)/z_sum;  %changed i to 1


%     %% Measure the expected errors
%     [~, D] = shortestpathtree(G, PATIENT);
%     overallcost = 0; 
%     for i = 1:1:NR_LOC
%         for j = 1:1:NR_LOC
%             approx_distance = sqrt((coordinate(j, 1) - coordinate(PATIENT, 1))^2 + (coordinate(j, 2) - coordinate(PATIENT, 2))^2 + (coordinate(j, 3) - coordinate(PATIENT, 3))^2); 
%             distance_error = abs(approx_distance - D(i)); 
%             costMatrix(i, j) = distance_error; 
%             overallcost = overallcost + distance_error * z(i, j);
%         end
%     end
%     overallcost = overallcost/NR_LOC; 
%     cost_distribution = cost_error_distribution(z, costMatrix, NR_LOC); 
end
