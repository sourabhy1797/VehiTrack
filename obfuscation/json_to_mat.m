% Read the JSON content from a file
jsonText = fileread('locationsets_new_100.json');

% Parse the JSON content into a MATLAB structure
jsonStruct = jsondecode(jsonText);

% Convert the structure to a cell structure using a recursive function
location_set = struct2cell_recursive(jsonStruct);

% Define a recursive function to convert the structure to a cell structure
function cellStruct = struct2cell_recursive(inputStruct)
    fieldNames = fieldnames(inputStruct);
    cellStruct = cell( numel(fieldNames), 1);

    for i = 1:numel(fieldNames)
        fieldName = fieldNames{i};
        fieldValue = inputStruct.(fieldName);

        if isstruct(fieldValue)
            cellStruct{i} = struct2cell_recursive(fieldValue);
        else
            cellStruct{i} = fieldValue;
        end
    end
end

% Now, cellStruct contains the JSON data in a cell structure
