function obfuscationMatrix = LPObfuscationMatrx(UL_matrix, distance_matrix, EPSILON)

%% Input
% UL_matrix: Utility loss matrix, each UL(i, k) represents the utility loss
% caused by the obfuscated location k given the real location i
% distance_matrix: each distance_matrix(i, j) represents the Haversine
% distance between location i and location j
% EPSILON: the privacy budget
%% Output
% the obfuscation matrix 


    NR_LOC = size(UL_matrix, 1); 

    %% Create the matrix for the Geo-indistinguishability constraints
    GeoI = sparse(NR_LOC*NR_LOC*(NR_LOC-1), NR_LOC*NR_LOC); 
    idx = 1; 
    for i = 1:1:NR_LOC
        for j = i+1:1:NR_LOC    
            for k = 1:1:NR_LOC
%                 GeoI(idx, (i-1)*NR_LOC + k) = exp(-EPSILON*distance_matrix(i, j)/2);
% 
%                 GeoI(idx, (j-1)*NR_LOC + k) = -exp(EPSILON*distance_matrix(i, j)/2);
% 
%                 idx = idx + 1;
% 
%                 GeoI(idx, (i-1)*NR_LOC + k) = -exp(EPSILON*distance_matrix(i, j)/2);
% 
%                 GeoI(idx, (j-1)*NR_LOC + k) = exp(-EPSILON*distance_matrix(i, j)/2);
% 
%                 idx = idx + 1;

                GeoI(idx, (i-1)*NR_LOC + k) = 1;
                GeoI(idx, (j-1)*NR_LOC + k) = -exp(EPSILON*distance_matrix(i, j));
                idx = idx + 1;
                GeoI(idx, (i-1)*NR_LOC + k) = -exp(EPSILON*distance_matrix(i, j));
                GeoI(idx, (j-1)*NR_LOC + k) = 1;
                idx = idx + 1;
            end
        end
    end
    b_GeoI = zeros(NR_LOC*NR_LOC*(NR_LOC-1), 1); 


    %% Create the cost vector for the objective function
    for i = 1:1:NR_LOC
        for k = 1:1:NR_LOC
            f((i-1)*NR_LOC + k) = UL_matrix(i, k); 
        end
    end

    %% Create the matrix for the probability unit measure constraints
    for i = 1:1:NR_LOC
        for j = 1:1:NR_LOC
            A_um(i, (i-1)*NR_LOC + j) = 1;
        end
    end  

    b_um = ones(NR_LOC, 1);


    %% Upper bound and lower bound of the decision variables
    lb = zeros(NR_LOC*NR_LOC, 1);
    ub = ones(NR_LOC*NR_LOC, 1);
    clear z;
    % z = {};

    [z, pfval, exitflag] = linprog(f, GeoI, b_GeoI, A_um, b_um, lb, ub);
    % [z{end+1}, pfval, exitflag] = linprog(f, GeoI, b_GeoI, A_um, b_um, lb, ub);
    if exitflag == -2
        options = optimoptions('linprog','Algorithm','interior-point');
        [z, pfval, exitflag] = linprog(f, GeoI, b_GeoI, A_um, b_um, lb, ub, options);
        % [z{end+1}, pfval, exitflag] = linprog(f, GeoI, b_GeoI, A_um, b_um, lb, ub);
    end


    obfuscationMatrix = reshape(z, NR_LOC, NR_LOC);
    obfuscationMatrix = obfuscationMatrix'; 

end