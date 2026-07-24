% Combine and Preprocess Multiple JSON Files

% List the names (or full paths) of the JSON files you want to combine.
jsonFiles = { 'OS_spatial_ev_master.json', ...
              'OS_spatial_ev_master2.json', ...
              'OS_spatial_ev_master3.json', ...
              'OS_spatial_ev_master3.1.json', ...
              'new.json'};

% Initialize an empty array to hold the combined data.
combinedData = [];

% Loop through each file and combine the data.
for k = 1:numel(jsonFiles)
    filename = jsonFiles{k};
    % Check if the file exists using 'exist'.
    if exist(filename, 'file') ~= 2
        error('File "%s" does not exist.', filename);
    end
    
    % Read the entire file content as text.
    jsonText = fileread(filename);
    
    % Decode the JSON text into a MATLAB variable.
    % Assuming each file is a JSON array of objects.
    data = jsondecode(jsonText);
    
    % Concatenate the arrays. (Ensure each file returns a compatible array.)
    combinedData = [combinedData; data];
end

% Preprocessing: Standardize the 'omega' field so that it is always a structure
% with fields 're' (real part) and 'im' (imaginary part). For entries where
% omega is a real number, set the imaginary part to 0.
for i = 1:numel(combinedData)
    if ~isstruct(combinedData(i).omega)
        % If omega is not a structure, assume it is a real number.
        omega_val = combinedData(i).omega;
        combinedData(i).omega = struct('re', omega_val, 'im', 0.0);
    else
        % If omega is a structure, ensure both 're' and 'im' fields exist.
        if ~isfield(combinedData(i).omega, 're')
            combinedData(i).omega.re = 0.0;
        end
        if ~isfield(combinedData(i).omega, 'im')
            combinedData(i).omega.im = 0.0;
        end
    end
end

% (Optional) Save the combined and standardized data into a new JSON file.
outputFilename = 'combined_data_standardized.json';
jsonCombinedText = jsonencode(combinedData);

% Open the file for writing.
fid = fopen(outputFilename, 'w');
if fid == -1
    error('Cannot open file "%s" for writing.', outputFilename);
end
fwrite(fid, jsonCombinedText, 'char');
fclose(fid);

fprintf('Combined and standardized JSON data saved to "%s".\n', outputFilename);

% Now "combinedData" contains the merged and preprocessed data.
% You can pass "combinedData" to your master analysis script.
