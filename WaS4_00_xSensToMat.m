%% XSENS DATA CONVERSION SCRIPT - WaS4 Study
% Converts Xsens .xlsx files to .mat format for faster loading
% Run this script after the pupilMerge preprocessing pipeline
%
% This script:
% 1. Finds all subjects in the data directory (matching pupilMerge structure)
% 2. Locates Xsens .xlsx files for each subject
% 3. Reads specified sheets from Excel files
% 4. Converts and saves data as .mat files in the raw data folder
% 5. Provides detailed progress reporting and error handling
%
% COMPATIBLE WITH: WaS4_00_pupilMerge.m pipeline

clear; clc; close all;

%% ==================== CONFIGURATION ====================

% Data paths - should match pupilMerge script
dataPath = '/Volumes/ergo/rohdaten/RohDaten_GRAIL/WaS4/DATA/';  % Main data directory (raw data)

% Processing options
processAllSubjects = true;  % Set to false to process specific subjects
specificSubjects = [];     % Only used if processAllSubjects = false
overwriteExisting = true;   % Set to true to reprocess existing .mat files

% Xsens file configuration
% Xsens file configuration
xsensSheets = {
    'General Information', 
    'Center of Mass', 
    'Segment Position', 
    'Joint Angles ZXY',
    'Segment Orientation - Quat',     % Quaternion data for gait detection
    'Segment Orientation - Euler'     % Euler angle data for gait detection
};
xsensVarNames = {
    'info', 
    'center_of_mass', 
    'segment_position', 
    'joint_angles',
    'segment_orientation_quat',       % Quaternion data for gait detection
    'segment_orientation_euler'       % Euler angle data for gait detection
};

% Output options
verboseOutput = true;
showProgress = true;

%% ==================== INITIALIZATION ====================

if verboseOutput
    fprintf('==========================================================\n');
    fprintf('XSENS DATA CONVERSION SCRIPT - WaS4 Study\n');
    fprintf('==========================================================\n\n');
    fprintf('📁 Data path: %s\n', dataPath);
    fprintf('📊 Sheets to convert: %s\n', strjoin(xsensSheets, ', '));
    fprintf('\n');
end

% Check if data directory exists
if ~exist(dataPath, 'dir')
    error('Data directory not found: %s', dataPath);
end

% Find all subject folders (matching pupilMerge pattern)
subjectFolders = dir(fullfile(dataPath, 'WaS_*'));
subjectFolders = subjectFolders([subjectFolders.isdir]);

if isempty(subjectFolders)
    error('No subject folders found in %s', dataPath);
end

% Filter subjects if not processing all
if ~processAllSubjects
    validFolders = {};
    for s = 1:length(specificSubjects)
        subjectNum = specificSubjects(s);
        folderPattern = sprintf('WaS_%03d', subjectNum);
        matchIdx = find(contains({subjectFolders.name}, folderPattern));
        if ~isempty(matchIdx)
            validFolders{end+1} = subjectFolders(matchIdx(1)).name;
        else
            warning('Subject %03d not found', subjectNum);
        end
    end
    subjectList = validFolders;
else
    subjectList = {subjectFolders.name};
end

if verboseOutput
    fprintf('🎯 Processing %d subjects:\n', length(subjectList));
    for i = 1:length(subjectList)
        fprintf('   - %s\n', subjectList{i});
    end
    fprintf('\n');
end

% Initialize tracking variables
processedSubjects = {};
failedSubjects = {};
skippedSubjects = {};
conversionStats = [];

% Define which sheets to process from XSens Excel files
sheetsToProcess = {
    'General Information',
    'Center of Mass', 
    'Segment Position',
    'Joint Angles ZXY',
    'Segment Orientation - Quat',
    'Segment Orientation - Euler'
};

% Initialize tracking variables
processedData = struct();
totalSheets = 0;
totalCells = 0;

%% ==================== MAIN PROCESSING LOOP ====================

for s = 1:length(subjectList)
    subjectFolder = subjectList{s};
    subjectPath = fullfile(dataPath, subjectFolder);
    
    % Extract subject number from folder name
    subjectNumStr = regexp(subjectFolder, 'WaS_(\d+)', 'tokens', 'once');
    if isempty(subjectNumStr)
        warning('Could not extract subject number from folder: %s', subjectFolder);
        failedSubjects{end+1} = sprintf('%s: Could not extract subject number', subjectFolder);
        continue;
    end
    subjectNum = str2double(subjectNumStr{1});
    
    if verboseOutput
        fprintf('==========================================================\n');
        fprintf('🔄 PROCESSING XSENS DATA: %s (Subject %d) [%d/%d]\n', subjectFolder, subjectNum, s, length(subjectList));
        fprintf('==========================================================\n\n');
    end
    
    try
        %% Check if already processed
        outputFile = fullfile(subjectPath, sprintf('%s_xsens.mat', subjectFolder));
        if exist(outputFile, 'file') && ~overwriteExisting
            if verboseOutput
                fprintf('⏭️  Already processed (use overwriteExisting=true to reprocess)\n');
                fprintf('   📄 Existing file: %s\n\n', outputFile);
            end
            skippedSubjects{end+1} = subjectFolder;
            continue;
        end
        
        %% Find Xsens Excel files
        if verboseOutput, fprintf('🔍 Searching for Xsens Excel files...\n'); end
        
        % Look for Excel files with subject number in filename
        xlsxPattern = sprintf('*%d*.xlsx', subjectNum);
        xlsxFiles = dir(fullfile(subjectPath, xlsxPattern));
        
        % Also try alternative patterns
        if isempty(xlsxFiles)
            xlsxPattern2 = sprintf('*%03d*.xlsx', subjectNum);
            xlsxFiles = dir(fullfile(subjectPath, xlsxPattern2));
        end
        
        % Generic Excel file search as fallback
        if isempty(xlsxFiles)
            xlsxFiles = dir(fullfile(subjectPath, '*.xlsx'));
        end
        
        if isempty(xlsxFiles)
            error('No Excel files found for subject %d in %s', subjectNum, subjectPath);
        end
        
        % Quick fix: Filter out system files
        xlsxFiles = xlsxFiles(~startsWith({xlsxFiles.name}, '._'));

        if verboseOutput
            fprintf('   📄 Found %d Excel file(s):\n', length(xlsxFiles));
            for f = 1:length(xlsxFiles)
                fileInfo = dir(fullfile(xlsxFiles(f).folder, xlsxFiles(f).name));
                fileSizeMB = fileInfo.bytes / (1024*1024);
                fprintf('      - %s (%.1f MB)\n', xlsxFiles(f).name, fileSizeMB);
            end
            fprintf('\n');
        end
        
        %% Initialize Xsens data structure
        if verboseOutput, fprintf('🔄 Converting Excel data to MATLAB format...\n'); end
        
        xsensData = struct();
        
        % Initialize all variable names
        for v = 1:length(xsensVarNames)
            xsensData.(xsensVarNames{v}) = {};
        end
        
        % Track conversion statistics
        totalSheetsProcessed = 0;
        totalCellsProcessed = 0;
        startTime = tic;
        
        %% Process each Excel file
        for f = 1:length(xlsxFiles)
            excelFile = fullfile(xlsxFiles(f).folder, xlsxFiles(f).name);
            
            if verboseOutput
                fprintf('   📊 Processing file %d/%d: %s\n', f, length(xlsxFiles), xlsxFiles(f).name);
            end
            
            % Check if file is accessible
            try
                % Try to get sheet names to verify file accessibility
                [~, sheetNames] = xlsfinfo(excelFile);
                if isempty(sheetNames)
                    warning('Could not read sheet names from file: %s', xlsxFiles(f).name);
                    continue;
                end
            catch ME
                warning('File access error for %s: %s', xlsxFiles(f).name, ME.message);
                continue;
            end
            
            %% Process each required sheet
            for sheetIdx = 1:length(sheetsToProcess)
                sheetName = sheetsToProcess{sheetIdx};
                
                if verboseOutput
                    fprintf('🔄 Reading sheet "%s"...\n', sheetName);
                end
                
                % Use the enhanced reading function
                data = read_large_excel_safe(excelFile, sheetName, verboseOutput);
                
                % FIXED: Store data in the correct xsensData structure
                if height(data) > 0
                    % Map sheet names to xsens variable names
                    switch lower(sheetName)
                        case 'general information'
                            varName = 'info';
                        case 'center of mass'
                            varName = 'center_of_mass';
                        case 'segment position'
                            varName = 'segment_position';
                        case 'joint angles zxy'
                            varName = 'joint_angles';
                        case 'segment orientation - quat'
                            varName = 'segment_orientation_quat';
                        case 'segment orientation - euler'
                            varName = 'segment_orientation_euler';
                        otherwise
                            varName = matlab.lang.makeValidName(lower(strrep(sheetName, ' ', '_')));
                    end
                    
                    % Store in xsensData structure (this gets saved!)
                    if ~isfield(xsensData, varName)
                        xsensData.(varName) = {};
                    end
                    xsensData.(varName){end+1} = data;
                    
                    % Update statistics
                    totalCellsProcessed = totalCellsProcessed + numel(data);
                    totalSheetsProcessed = totalSheetsProcessed + 1;
                    
                    if verboseOutput
                        fprintf('   ✅ Stored %d rows × %d cols in xsensData.%s\n', ...
                               height(data), width(data), varName);
                    end
                else
                    if verboseOutput
                        fprintf('   ⚠️ No data in sheet "%s"\n', sheetName);
                    end
                end
            end
        end
        
        processingTime = toc(startTime);
        
        %% Save converted data
        if verboseOutput, fprintf('\n💾 Saving converted data...\n'); end
        
        % Add metadata to the saved structure
        conversionInfo = struct();
        conversionInfo.conversion_date = datetime('now');
        conversionInfo.original_files = {xlsxFiles.name};
        conversionInfo.sheets_converted = xsensSheets;
        conversionInfo.processing_time_seconds = processingTime;
        conversionInfo.total_sheets_processed = totalSheetsProcessed;
        conversionInfo.total_cells_processed = totalCellsProcessed;
        conversionInfo.matlab_version = version;
        conversionInfo.script_name = mfilename;
        
        % Validate required data for gait detection
        requiredFields = {'segment_orientation_quat', 'segment_orientation_euler'};
        missingFields = [];

        for f = 1:length(requiredFields)
            if ~isfield(xsensData, requiredFields{f}) || isempty(xsensData.(requiredFields{f}))
                missingFields{end+1} = requiredFields{f};
            end
        end

        if ~isempty(missingFields)
            warning('Missing required data for gait detection: %s', strjoin(missingFields, ', '));
            
            % Add empty placeholder if data is missing
            for f = 1:length(missingFields)
                if ~isfield(xsensData, missingFields{f})
                    xsensData.(missingFields{f}) = {};
                end
            end
        end

        % Add sampling rate information if available
        if ~isempty(xsensData.info)
            try
                infoTable = xsensData.info{1};
                frameRateRow = contains(infoTable.Parameter, 'Frame Rate', 'IgnoreCase', true);
                if any(frameRateRow)
                    xsensData.sampling_rate = str2double(infoTable.Value{frameRateRow});
                end
            catch
                warning('Could not extract sampling rate from info sheet');
            end
        end

        % Save with compression for large files
        save(outputFile, 'xsensData', 'conversionInfo', '-v7.3');
        
        % Verify saved file
        savedFileInfo = dir(outputFile);
        savedFileSizeMB = savedFileInfo.bytes / (1024*1024);
        
        if verboseOutput
            fprintf('   ✅ Data saved: %s\n', outputFile);
            fprintf('   📊 File size: %.1f MB\n', savedFileSizeMB);
            fprintf('   ⏱️  Processing time: %.1f seconds\n', processingTime);
            fprintf('   📈 Performance: %d sheets, %d cells processed\n', ...
                   totalSheetsProcessed, totalCellsProcessed);
        end
        
        % Store statistics
        stats = struct();
        stats.subject = subjectFolder;
        stats.excel_files_found = length(xlsxFiles);
        stats.sheets_processed = totalSheetsProcessed;
        stats.cells_processed = totalCellsProcessed;
        stats.processing_time = processingTime;
        stats.output_file_size_mb = savedFileSizeMB;
        conversionStats = [conversionStats; stats];
        
        % Mark as successfully processed
        processedSubjects{end+1} = subjectFolder;
        
        if verboseOutput
            fprintf('\n✅ %s CONVERSION COMPLETED!\n', subjectFolder);
            fprintf('   📁 Excel files processed: %d\n', length(xlsxFiles));
            fprintf('   📊 Total sheets: %d, Total cells: %d\n', totalSheetsProcessed, totalCellsProcessed);
            fprintf('   ⏱️  Processing time: %.1f seconds\n', processingTime);
        end
        
    catch ME
        % Handle errors gracefully
        failedSubjects{end+1} = sprintf('%s: %s', subjectFolder, ME.message);
        
        if verboseOutput
            fprintf('\n❌ ERROR processing %s:\n', subjectFolder);
            fprintf('   %s\n', ME.message);
            if length(ME.stack) > 0
                fprintf('   Line %d in %s\n', ME.stack(1).line, ME.stack(1).name);
            end
        end
        
        % Continue with next subject
        continue;
    end
end

%% ==================== FINAL SUMMARY ====================

if verboseOutput
    fprintf('\n==========================================================\n');
    fprintf('XSENS CONVERSION SUMMARY\n');
    fprintf('==========================================================\n\n');
    
    fprintf('✅ SUCCESSFULLY PROCESSED (%d subjects):\n', length(processedSubjects));
    for i = 1:length(processedSubjects)
        if i <= length(conversionStats)
            fprintf('   %s (%d files, %d sheets, %.1f MB, %.1fs)\n', processedSubjects{i}, ...
                   conversionStats(i).excel_files_found, conversionStats(i).sheets_processed, ...
                   conversionStats(i).output_file_size_mb, conversionStats(i).processing_time);
        else
            fprintf('   %s\n', processedSubjects{i});
        end
    end
    
    if ~isempty(skippedSubjects)
        fprintf('\n⏭️  SKIPPED (already processed) (%d subjects):\n', length(skippedSubjects));
        for i = 1:length(skippedSubjects)
            fprintf('   %s\n', skippedSubjects{i});
        end
    end
    
    if ~isempty(failedSubjects)
        fprintf('\n❌ FAILED (%d subjects):\n', length(failedSubjects));
        for i = 1:length(failedSubjects)
            fprintf('   %s\n', failedSubjects{i});
        end
    end
    
    % Overall statistics
    if ~isempty(conversionStats)
        totalFiles = sum([conversionStats.excel_files_found]);
        totalSheets = sum([conversionStats.sheets_processed]);
        totalCells = sum([conversionStats.cells_processed]);
        totalTime = sum([conversionStats.processing_time]);
        totalSize = sum([conversionStats.output_file_size_mb]);
        
        fprintf('\n📊 CONVERSION STATISTICS:\n');
        fprintf('   Excel files processed: %d\n', totalFiles);
        fprintf('   Sheets converted: %d\n', totalSheets);
        fprintf('   Cells processed: %,d\n', totalCells);
        fprintf('   Total processing time: %.1f seconds (%.1f minutes)\n', totalTime, totalTime/60);
        fprintf('   Total output size: %.1f MB\n', totalSize);
        fprintf('   Average processing rate: %,.0f cells/second\n', totalCells/totalTime);
        
        if length(processedSubjects) > 1
            fprintf('\n📈 PER-SUBJECT AVERAGES:\n');
            fprintf('   Files per subject: %.1f\n', totalFiles / length(processedSubjects));
            fprintf('   Sheets per subject: %.1f\n', totalSheets / length(processedSubjects));
            fprintf('   Processing time per subject: %.1f seconds\n', totalTime / length(processedSubjects));
            fprintf('   Output size per subject: %.1f MB\n', totalSize / length(processedSubjects));
        end
    end
    
    fprintf('\n💾 Data directory: %s\n', dataPath);
    fprintf('📋 Conversion completed: %s\n', datestr(now));
    fprintf('==========================================================\n');
end

%% ==================== HELPER FUNCTIONS ====================

%% Replace the helper function at the bottom of your script with this:

%% Replace the helper function at the bottom of your script with this:

function data = read_large_excel_safe(filepath, sheetName, verboseOutput)
    % Enhanced Excel reader for large files with local copy strategy
    % Solution 1: Copy to local drive first for better performance
    
    data = table();
    
    % Check if file exists
    if ~exist(filepath, 'file')
        if verboseOutput
            fprintf('     ❌ File not found: %s\n', filepath);
        end
        return;
    end
    
    % Get file size for strategy decision
    try
        fileInfo = dir(filepath);
        fileSizeMB = fileInfo.bytes / (1024^2);
    catch
        if verboseOutput
            fprintf('     ⚠️ Could not determine file size\n');
        end
        fileSizeMB = 0;
    end
    
    if verboseOutput
        fprintf('     📊 File size: %.1f MB\n', fileSizeMB);
    end
    
    % Strategy 1: For large files (>50MB) or network drives, copy locally first
    isNetworkDrive = contains(filepath, '/Volumes/') || contains(filepath, '\\') || contains(filepath, 'smb://');
    useLocalCopy = fileSizeMB > 50 || isNetworkDrive;
    
    if useLocalCopy && verboseOutput
        fprintf('     📋 Using local copy strategy (Large file or network drive)\n');
    end
    
    originalFile = filepath;
    cleanupRequired = false;
    
    try
        if useLocalCopy
            % Create temporary local copy
            [~, filename, ext] = fileparts(filepath);
            tempDir = tempdir;
            localFile = fullfile(tempDir, [filename '_temp' ext]);
            
            if verboseOutput
                fprintf('     🔄 Copying to local temp: %s\n', localFile);
            end
            
            % Copy file with progress indication for very large files
            if fileSizeMB > 500
                if verboseOutput
                    fprintf('     ⏳ Large file copy in progress (this may take a while)...\n');
                end
            end
            
            copyfile(originalFile, localFile);
            filepath = localFile;  % Use local copy
            cleanupRequired = true;
            
            if verboseOutput
                fprintf('     ✅ File copied to local drive\n');
            end
        end
        
        % Strategy 2: Try multiple reading methods in order of preference
        methods = {'readtable', 'readmatrix', 'xlsread_fallback'};
        
        for methodIdx = 1:length(methods)
            method = methods{methodIdx};
            
            try
                if verboseOutput && methodIdx > 1
                    fprintf('     🔄 Trying method %d: %s\n', methodIdx, method);
                end
                
                switch method
                    case 'readtable'
                        % Standard method
                        data = readtable(filepath, 'Sheet', sheetName, 'ReadVariableNames', true);
                        
                    case 'readmatrix'
                        % Fallback for problematic files
                        if verboseOutput
                            fprintf('     🔧 Using readmatrix fallback\n');
                        end
                        rawData = readmatrix(filepath, 'Sheet', sheetName);
                        
                        if ~isempty(rawData)
                            % Create appropriate column names based on sheet
                            if strcmpi(sheetName, 'Center of Mass')
                                colNames = {'Frame', 'CoM_pos_x', 'CoM_pos_y', 'CoM_pos_z', 'CoM_vel_x', 'CoM_vel_y', 'CoM_vel_z', ...
                                           'CoM_acc_x', 'CoM_acc_y', 'CoM_acc_z'};
                            elseif strcmpi(sheetName, 'General Information')
                                colNames = {'Parameter', 'Value'};
                            elseif strcmpi(sheetName, 'Segment Position')
                                colNames = {'Frame', 'Time'};
                                % Add segment position columns dynamically
                                for i = 3:size(rawData, 2)
                                    colNames{end+1} = sprintf('Segment%d', i-2);
                                end
                            elseif strcmpi(sheetName, 'Joint Angles ZXY')
                                colNames = {'Frame', 'Time'};
                                % Add joint angle columns dynamically
                                for i = 3:size(rawData, 2)
                                    colNames{end+1} = sprintf('Joint%d', i-2);
                                end
                            else
                                % Generic column names
                                colNames = arrayfun(@(x) sprintf('Var%d', x), 1:size(rawData,2), 'UniformOutput', false);
                            end
                            
                            % Ensure we don't have more columns than data
                            numCols = min(length(colNames), size(rawData,2));
                            colNames = colNames(1:numCols);
                            
                            % Convert to table
                            data = array2table(rawData(:, 1:numCols), 'VariableNames', colNames);
                        end
                        
                    case 'xlsread_fallback'
                        % Last resort: older xlsread function
                        if verboseOutput
                            fprintf('     🔧 Using xlsread fallback (last resort)\n');
                        end
                        
                        try
                            [num, txt, raw] = xlsread(filepath, sheetName);
                            
                            if ~isempty(num)
                                % Get proper column names based on sheet type
                                colNames = get_xsens_column_names(sheetName, size(num, 2));
                                
                                % If we have text data, try to extract headers
                                if ~isempty(txt) && size(txt, 1) > 0
                                    % Look for header row in text data
                                    headerRow = [];
                                    for row = 1:min(3, size(txt, 1)) % Check first 3 rows for headers
                                        if sum(~cellfun(@isempty, txt(row, :))) > size(num, 2) * 0.5
                                            headerRow = txt(row, :);
                                            break;
                                        end
                                    end
                                    
                                    if ~isempty(headerRow)
                                        % Clean and validate headers
                                        validHeaders = {};
                                        for h = 1:min(length(headerRow), size(num, 2))
                                            if ~isempty(headerRow{h}) && (ischar(headerRow{h}) || isstring(headerRow{h}))
                                                cleanHeader = strtrim(char(headerRow{h}));
                                                cleanHeader = regexprep(cleanHeader, '[^\w]', '_'); % Replace special chars
                                                cleanHeader = matlab.lang.makeValidName(cleanHeader);
                                                if ~isempty(cleanHeader) && ~strcmp(cleanHeader, 'x_')
                                                    validHeaders{h} = cleanHeader;
                                                else
                                                    validHeaders{h} = colNames{h};
                                                end
                                            else
                                                validHeaders{h} = colNames{h};
                                            end
                                        end
                                        
                                        % Use extracted headers if we have enough valid ones
                                        if length(validHeaders) >= size(num, 2) * 0.7
                                            colNames = validHeaders(1:size(num, 2));
                                            if verboseOutput
                                                fprintf('     📋 Extracted headers from Excel file\n');
                                            end
                                        end
                                    end
                                end
                                
                                % Create table with proper column names
                                data = array2table(num, 'VariableNames', colNames(1:size(num, 2)));
                                
                            elseif ~isempty(raw)
                                % Use raw data and try to convert
                                if size(raw, 1) > 1
                                    % Try to use first row as headers
                                    headers = raw(1, :);
                                    dataRows = raw(2:end, :);
                                    
                                    % Convert headers to valid variable names
                                    validHeaders = {};
                                    hasValidHeaders = false;
                                    
                                    for h = 1:length(headers)
                                        if ~isempty(headers{h}) && (ischar(headers{h}) || isstring(headers{h}))
                                            cleanHeader = strtrim(char(headers{h}));
                                            if ~isnumeric(headers{h}) && length(cleanHeader) > 0
                                                cleanHeader = regexprep(cleanHeader, '[^\w]', '_');
                                                validHeaders{h} = matlab.lang.makeValidName(cleanHeader);
                                                hasValidHeaders = true;
                                            else
                                                validHeaders{h} = sprintf('Var%d', h);
                                            end
                                        else
                                            validHeaders{h} = sprintf('Var%d', h);
                                        end
                                    end
                                    
                                    % If no valid headers found, use sheet-specific names
                                    if ~hasValidHeaders
                                        validHeaders = get_xsens_column_names(sheetName, length(validHeaders));
                                    end
                                    
                                    % Try to convert data to numeric
                                    try
                                        numericData = cell2mat(dataRows);
                                        if ~isempty(numericData)
                                            data = array2table(numericData, 'VariableNames', validHeaders(1:size(numericData,2)));
                                        end
                                    catch
                                        % Mixed data types, create table directly from cell array
                                        if size(dataRows, 2) <= length(validHeaders)
                                            data = cell2table(dataRows, 'VariableNames', validHeaders(1:size(dataRows,2)));
                                        end
                                    end
                                end
                            end
                        catch
                            % xlsread also failed
                            data = table();
                        end
                end
                
                % Check if we successfully read data
                if height(data) > 0
                    if verboseOutput
                        fprintf('     ✅ Successfully read %d rows, %d columns using %s\n', ...
                               height(data), width(data), method);
                    end
                    break;  % Success, exit method loop
                else
                    if verboseOutput && methodIdx < length(methods)
                        fprintf('     ⚠️ No data returned, trying next method\n');
                    end
                end
                
            catch ME
                if verboseOutput
                    if methodIdx < length(methods)
                        fprintf('     ⚠️ Method %s failed: %s\n', method, ME.message);
                        fprintf('     🔄 Trying next method...\n');
                    else
                        fprintf('     ❌ All methods failed. Last error: %s\n', ME.message);
                    end
                end
                
                if methodIdx == length(methods)
                    % All methods failed
                    data = table();
                end
            end
        end
        
    catch ME
        if verboseOutput
            fprintf('     ❌ Critical error in file processing: %s\n', ME.message);
        end
        data = table();
    end
    
    % Cleanup temporary file
    if cleanupRequired && exist(filepath, 'file')
        try
            delete(filepath);
            if verboseOutput
                fprintf('     🗑️ Cleaned up temporary file\n');
            end
        catch
            if verboseOutput
                fprintf('     ⚠️ Could not clean up temporary file: %s\n', filepath);
            end
        end
    end
    
    % Final validation
    if height(data) == 0
        if verboseOutput
            fprintf('     ❌ No data could be read from sheet "%s"\n', sheetName);
        end
    end
end

function colNames = get_xsens_column_names(sheetName, numCols)
    % Get appropriate column names for XSens sheets based on sheet type
    
    switch lower(strtrim(sheetName))
        case 'center of mass'
            colNames = {
                'Frame', 'CoM_pos_x', 'CoM_pos_y', 'CoM_pos_z', ...
                'CoM_vel_x', 'CoM_vel_y', 'CoM_vel_z', ...
                'CoM_acc_x', 'CoM_acc_y', 'CoM_acc_z'
            };
            
        case 'general information'
            colNames = {'Parameter', 'Value'};
            
        case 'segment position'
            colNames = {'Frame', 'Time'};
            % Add segment-specific columns
            segmentNames = {
                'Pelvis', 'L5', 'L3', 'T12', 'T8', 'Neck', 'Head', ...
                'RightShoulder', 'RightUpperArm', 'RightForearm', 'RightHand', ...
                'LeftShoulder', 'LeftUpperArm', 'LeftForearm', 'LeftHand', ...
                'RightUpperLeg', 'RightLowerLeg', 'RightFoot', 'RightToe', ...
                'LeftUpperLeg', 'LeftLowerLeg', 'LeftFoot', 'LeftToe'
            };
            
            for seg = 1:length(segmentNames)
                if length(colNames) >= numCols, break; end
                for coord = {'X', 'Y', 'Z'}
                    colNames{end+1} = sprintf('%s_%s', segmentNames{seg}, coord{1});
                    if length(colNames) >= numCols, break; end
                end
            end
            
        case 'joint angles zxy'
            colNames = {'Frame', 'Time'};
            % Add joint-specific columns
            jointNames = {
                'jL5S1', 'jL4L3', 'jL1T12', 'jT9T8', 'jT1C7', 'jC1Head', ...
                'jRightT4Shoulder', 'jRightShoulder', 'jRightElbow', 'jRightWrist', ...
                'jLeftT4Shoulder', 'jLeftShoulder', 'jLeftElbow', 'jLeftWrist', ...
                'jRightHip', 'jRightKnee', 'jRightAnkle', 'jRightBallFoot', ...
                'jLeftHip', 'jLeftKnee', 'jLeftAnkle', 'jLeftBallFoot'
            };
            
            for joint = 1:length(jointNames)
                if length(colNames) >= numCols, break; end
                for angle = {'Z', 'X', 'Y'}
                    colNames{end+1} = sprintf('%s_%s', jointNames{joint}, angle{1});
                    if length(colNames) >= numCols, break; end
                end
            end

        case 'angular kinematics'
            colNames = {'Frame', 'Time'};
            segmentNames = {
                'Pelvis', 'L5', 'L3', 'T12', 'T8', 'Neck', 'Head', ...
                'RightShoulder', 'RightUpperArm', 'RightForearm', 'RightHand', ...
                'LeftShoulder', 'LeftUpperArm', 'LeftForearm', 'LeftHand', ...
                'RightUpperLeg', 'RightLowerLeg', 'RightFoot', 'RightToe', ...
                'LeftUpperLeg', 'LeftLowerLeg', 'LeftFoot', 'LeftToe'
            };
            
            for seg = 1:length(segmentNames)
                if length(colNames) >= numCols, break; end
                for coord = {'AngVel_X', 'AngVel_Y', 'AngVel_Z', 'AngAcc_X', 'AngAcc_Y', 'AngAcc_Z'}
                    colNames{end+1} = sprintf('%s_%s', segmentNames{seg}, coord{1});
                    if length(colNames) >= numCols, break; end
                end
            end
            
        case 'segment kinematics'
            colNames = {'Frame', 'Time'};
            segmentNames = {
                'Pelvis', 'L5', 'L3', 'T12', 'T8', 'Neck', 'Head', ...
                'RightShoulder', 'RightUpperArm', 'RightForearm', 'RightHand', ...
                'LeftShoulder', 'LeftUpperArm', 'LeftForearm', 'LeftHand', ...
                'RightUpperLeg', 'RightLowerLeg', 'RightFoot', 'RightToe', ...
                'LeftUpperLeg', 'LeftLowerLeg', 'LeftFoot', 'LeftToe'
            };
            
            for seg = 1:length(segmentNames)
                if length(colNames) >= numCols, break; end
                for metric = {'Acc_X', 'Acc_Y', 'Acc_Z', 'Vel_X', 'Vel_Y', 'Vel_Z'}
                    colNames{end+1} = sprintf('%s_%s', segmentNames{seg}, metric{1});
                    if length(colNames) >= numCols, break; end
                end
            end

        case 'segment orientation - quat'
            colNames = {'Frame', 'Time'};
            segmentNames = {
                'Pelvis', 'L5', 'L3', 'T12', 'T8', 'Neck', 'Head', ...
                'RightShoulder', 'RightUpperArm', 'RightForearm', 'RightHand', ...
                'LeftShoulder', 'LeftUpperArm', 'LeftForearm', 'LeftHand', ...
                'RightUpperLeg', 'RightLowerLeg', 'RightFoot', 'RightToe', ...
                'LeftUpperLeg', 'LeftLowerLeg', 'LeftFoot', 'LeftToe'
            };
            
            for seg = 1:length(segmentNames)
                if length(colNames) >= numCols, break; end
                for quat = {'q0', 'q1', 'q2', 'q3'}
                    colNames{end+1} = sprintf('%s_%s', segmentNames{seg}, quat{1});
                    if length(colNames) >= numCols, break; end
                end
            end

        case 'segment orientation - euler'
            colNames = {'Frame', 'Time'};
            segmentNames = {
                'Pelvis', 'L5', 'L3', 'T12', 'T8', 'Neck', 'Head', ...
                'RightShoulder', 'RightUpperArm', 'RightForearm', 'RightHand', ...
                'LeftShoulder', 'LeftUpperArm', 'LeftForearm', 'LeftHand', ...
                'RightUpperLeg', 'RightLowerLeg', 'RightFoot', 'RightToe', ...
                'LeftUpperLeg', 'LeftLowerLeg', 'LeftFoot', 'LeftToe'
            };
            
            for seg = 1:length(segmentNames)
                if length(colNames) >= numCols, break; end
                for angle = {'X', 'Y', 'Z'}
                    colNames{end+1} = sprintf('%s_%s', segmentNames{seg}, angle{1});
                    if length(colNames) >= numCols, break; end
                end
            end

        otherwise
            % Generic column names as fallback
            colNames = arrayfun(@(x) sprintf('Var%d', x), 1:numCols, 'UniformOutput', false);

    end
    
    % Ensure we have exactly the right number of columns
    if length(colNames) < numCols
        % Add more generic columns if needed
        for i = (length(colNames) + 1):numCols
            colNames{i} = sprintf('Col%d', i);
        end
    elseif length(colNames) > numCols
        % Trim if we have too many
        colNames = colNames(1:numCols);
    end
    
    % Ensure all names are valid MATLAB variable names
    colNames = cellfun(@(x) matlab.lang.makeValidName(x), colNames, 'UniformOutput', false);
end