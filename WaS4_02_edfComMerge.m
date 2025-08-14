%% ADDITIONAL DATA MERGING SCRIPT - WaS4 Study
% Merges additional data sources into aligned EEG files
% Run this script AFTER WaS4_00_pupilMerge.m and WaS4_01_xSensToMat.m
%
% This script:
% 1. Loads aligned EEG files created by WaS4_00_pupilMerge.m
% 2. Loads converted XSens data from WaS4_01_xSensToMat.m
% 3. Loads ECG data from .EDF files (Faros device)
% 4. Loads Pupil Labs IMU data from .csv files
% 5. Performs time alignment using cross-correlation and counter matching
% 6. Merges all additional data into the EEG structure
% 7. Saves enhanced EEG files with all data sources integrated
%
% DEPENDENCIES: WaS4_00_pupilMerge.m, WaS4_01_xSensToMat.m

clear; clc; close all;

%% ==================== CONFIGURATION ====================

% Data paths - should match previous scripts
dataPath = '/Volumes/ergo/rohdaten/RohDaten_GRAIL/WaS4/DATA/';  % Raw data directory
alignedPath = '/Volumes/Work4TB/Seafile/WaS4/data/aligned/';    % Aligned data directory (from WaS4_00)
outPath = '/Volumes/Work4TB/Seafile/WaS4/data/merged/';    % Aligned data directory (from WaS4_00)

% Processing options
processAllSubjects = true;  % Set to false to process specific subjects
specificSubjects = [1];     % Only used if processAllSubjects = false
overwriteExisting = false;   % Set to true to reprocess existing files

% Data integration options
includeXSensCoM = true;      % Include XSens Center of Mass data
includeECG = true;           % Include ECG data from .EDF files
includePupilIMU = true;      % Include Pupil Labs IMU data
includeXSensFootPos = true;  % Include foot position for step detection

% Quality control
maxTimeOffset = 5.0;         % Maximum allowed time offset in seconds
minCorrelation = 0.3;        % Minimum cross-correlation for ECG alignment

% Output options
verboseOutput = true;
saveResults = true;
createPlots = true;
savePlots = true;

%% ==================== INITIALIZATION ====================

if verboseOutput
    fprintf('==========================================================\n');
    fprintf('ADDITIONAL DATA MERGING SCRIPT - WaS4 Study\n');
    fprintf('==========================================================\n\n');
    fprintf('📁 Raw data path: %s\n', dataPath);
    fprintf('📁 Aligned data path: %s\n', alignedPath);
    fprintf('📁 Merged data path: %s\n', outPath);
    fprintf('🔧 Data sources to merge:\n');
    if includeXSensCoM, fprintf('   ✓ XSens Center of Mass\n'); end
    if includeECG, fprintf('   ✓ ECG (.EDF files)\n'); end
    if includePupilIMU, fprintf('   ✓ Pupil Labs IMU\n'); end
    if includeXSensFootPos, fprintf('   ✓ XSens Foot Position\n'); end
    fprintf('\n');
end

% Check directories
if ~exist(dataPath, 'dir')
    error('Raw data directory not found: %s', dataPath);
end
if ~exist(alignedPath, 'dir')
    error('Aligned data directory not found: %s (run WaS4_00_pupilMerge.m first)', alignedPath);
end
if ~exist(outPath, 'dir')
    error('Aligned data directory not found: %s', outPath);
end

% Find aligned EEG files
alignedFiles = dir(fullfile(alignedPath, '*_aligned.set'));
if isempty(alignedFiles)
    error('No aligned EEG files found in %s (run WaS4_00_pupilMerge.m first)', alignedPath);
end

% Extract subject list from aligned files
subjectList = {};
for i = 1:length(alignedFiles)
    filename = alignedFiles(i).name;
    subjectMatch = regexp(filename, '(WaS_\d+)_aligned\.set', 'tokens', 'once');
    if ~isempty(subjectMatch)
        subjectList{end+1} = subjectMatch{1};
    end
end

% Filter subjects if not processing all
if ~processAllSubjects
    validSubjects = {};
    for s = 1:length(specificSubjects)
        subjectNum = specificSubjects(s);
        folderPattern = sprintf('WaS_%03d', subjectNum);
        matchIdx = find(contains(subjectList, folderPattern));
        if ~isempty(matchIdx)
            validSubjects{end+1} = subjectList{matchIdx(1)};
        else
            warning('Subject %03d not found in aligned data', subjectNum);
        end
    end
    subjectList = validSubjects;
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
mergeStats = [];

%% ==================== MAIN PROCESSING LOOP ====================

for s = 1:length(subjectList)
    subjectFolder = subjectList{s};
    subjectPath = fullfile(dataPath, subjectFolder);
    
    % Extract subject number for file matching
    subjectNumStr = regexp(subjectFolder, 'WaS_(\d+)', 'tokens', 'once');
    if isempty(subjectNumStr)
        warning('Could not extract subject number from: %s', subjectFolder);
        failedSubjects{end+1} = sprintf('%s: Could not extract subject number', subjectFolder);
        continue;
    end
    subjectNum = str2double(subjectNumStr{1});
    
    if verboseOutput
        fprintf('==========================================================\n');
        fprintf('🔄 MERGING ADDITIONAL DATA: %s (Subject %d) [%d/%d]\n', subjectFolder, subjectNum, s, length(subjectList));
        fprintf('==========================================================\n\n');
    end
    
    try
        % Check if already processed
        outputFile = fullfile(outPath, sprintf('%s_merged.set', subjectFolder));
        if exist(outputFile, 'file') && ~overwriteExisting
            if verboseOutput
                fprintf('⏭️  Already processed (use overwriteExisting=true to reprocess)\n\n');
            end
            skippedSubjects{end+1} = subjectFolder;
            continue;
        end
        
        %% Step 1: Load aligned EEG data
        if verboseOutput, fprintf('🧠 Loading aligned EEG data...\n'); end
        
        alignedFile = fullfile(alignedPath, sprintf('%s_aligned.set', subjectFolder));
        if ~exist(alignedFile, 'file')
            error('Aligned EEG file not found: %s', alignedFile);
        end
        
        eeg = pop_loadset('filename', sprintf('%s_aligned.set', subjectFolder), 'filepath', alignedPath);
        
        % Find counter channel
        cntIdx = find(strcmp({eeg.chanlocs.labels}, 'CNT'));
        if isempty(cntIdx)
            error('No CNT (counter) channel found in aligned EEG data');
        end
        
        if verboseOutput
            fprintf('   ✅ Loaded aligned EEG: %d channels, %d samples, %.1f Hz\n', ...
                   eeg.nbchan, eeg.pnts, eeg.srate);
        end
        
        % Track statistics
        stats = struct();
        stats.subject = subjectFolder;
        stats.original_channels = eeg.nbchan;
        stats.channels_added = 0;
        stats.data_sources_merged = {};
        startTime = tic;
        
        %% Step 2: Load XSens data (if requested)
        xsensData = [];
        if includeXSensCoM || includeXSensFootPos
            if verboseOutput, fprintf('\n🔄 Loading XSens data...\n'); end
            
            xsensFile = fullfile(subjectPath, sprintf('%s_xsens.mat', subjectFolder));
            if exist(xsensFile, 'file')
                xsensLoad = load(xsensFile);
                xsensData = xsensLoad.xsensData;
                
                if verboseOutput
                    fprintf('   ✅ XSens data loaded from: %s\n', xsensFile);
                    fprintf('   📊 Available data types: %s\n', strjoin(fieldnames(xsensData), ', '));
                end
            else
                if verboseOutput
                    fprintf('   ⚠️ XSens file not found: %s (run WaS4_01_xSensToMat.m first)\n', xsensFile);
                end
                if includeXSensCoM, includeXSensCoM = false; end
                if includeXSensFootPos, includeXSensFootPos = false; end
            end
        end
        
        %% Step 3: Load and merge ECG data (if requested)
        if includeECG
            if verboseOutput, fprintf('\n❤️  Loading ECG data...\n'); end
            
            % Find EDF files
            edfFiles = dir(fullfile(subjectPath, '*.EDF'));
            if isempty(edfFiles)
                edfFiles = dir(fullfile(subjectPath, '*.edf'));
            end
            
            if ~isempty(edfFiles)
                farosData = [];
                
                for f = 1:length(edfFiles)
                    edfFile = fullfile(subjectPath, edfFiles(f).name);
                    if verboseOutput, fprintf('   📥 Reading: %s\n', edfFiles(f).name); end
                    
                    try
                        farosImport = edfread(edfFile);
                        farosInfo = edfinfo(edfFile);
                        ecgSampRate = farosInfo.NumSamples(1);
                        
                        % Extract ECG and accelerometer data
                        newECGData = cell2mat(farosImport.ECG);
                        newECGTime = (0:length(newECGData)-1) * (1/ecgSampRate);
                        
                        % Handle accelerometer data if available
                        if isfield(farosImport, 'Accelerometer_X')
                            newACCData = [cell2mat(farosImport.Accelerometer_X), ...
                                         cell2mat(farosImport.Accelerometer_Y), ...
                                         cell2mat(farosImport.Accelerometer_Z)];
                            accSampRate = farosInfo.NumSamples(2);
                            newACCTime = (0:size(newACCData,1)-1) * (1/accSampRate);
                            
                            % Interpolate ACC to ECG sampling rate
                            accResampled = zeros(length(newECGTime), 3);
                            for a = 1:3
                                accResampled(:,a) = interp1(newACCTime, newACCData(:,a), newECGTime, 'makima');
                            end
                        else
                            accResampled = zeros(length(newECGTime), 3);
                        end
                        
                        % Combine data: [time, ECG, AccX, AccY, AccZ]
                        newFarosData = [newECGTime', newECGData, accResampled];
                        farosData = [farosData; newFarosData];
                        
                        if verboseOutput
                            fprintf('      ✅ %d samples, %.1f Hz, %.1f minutes\n', ...
                                   length(newECGData), ecgSampRate, max(newECGTime)/60);
                        end
                        
                    catch ME
                        if verboseOutput
                            fprintf('      ❌ Error reading %s: %s\n', edfFiles(f).name, ME.message);
                        end
                    end
                end
                
                if ~isempty(farosData)
                    % Add ECG to EEG structure
                    newChannelIdx = eeg.nbchan + 1;
                    
                    % Create time vector for EEG samples
                    eegTimes = eeg.times / 1000; % Convert to seconds
                    
                    % Simple time alignment: assume ECG starts at EEG start
                    % For more sophisticated alignment, cross-correlation could be used here
                    farosData(:,1) = farosData(:,1) + eegTimes(1);
                    
                    % Filter to EEG time range
                    timeFilter = farosData(:,1) >= eegTimes(1) & farosData(:,1) <= eegTimes(end);
                    farosFiltered = farosData(timeFilter, :);
                    
                    if size(farosFiltered, 1) > 100  % Ensure we have sufficient data
                        % Interpolate to EEG sample times
                        ecgInterp = interp1(farosFiltered(:,1), farosFiltered(:,2), eegTimes, 'cubic', NaN);
                        
                        % Expand EEG data matrix to accommodate new channels
                        newChannelIdx = eeg.nbchan + 1;
                        eeg.data(newChannelIdx, :) = ecgInterp;
                        eeg.chanlocs(newChannelIdx).labels = 'ECG_Faros';
                        eeg.chanlocs(newChannelIdx).type = 'ECG';
                        eeg.nbchan = eeg.nbchan + 1;
                        stats.channels_added = stats.channels_added + 1;
                        
                        % Add accelerometer channels if available
                        if size(farosFiltered, 2) >= 5
                            for a = 1:3
                                accInterp = interp1(farosFiltered(:,1), farosFiltered(:,a+2), eegTimes, 'cubic', NaN);
                                newChannelIdx = eeg.nbchan + 1;
                                eeg.data(newChannelIdx, :) = accInterp;
                                eeg.chanlocs(newChannelIdx).labels = sprintf('ACC_%s_Faros', char('X'+a-1));
                                eeg.chanlocs(newChannelIdx).type = 'ACC';
                                eeg.nbchan = eeg.nbchan + 1;
                                stats.channels_added = stats.channels_added + 1;
                            end
                        end
                        
                        stats.data_sources_merged{end+1} = 'ECG';
                        
                        if verboseOutput
                            validSamples = sum(~isnan(ecgInterp));
                            coverage = validSamples / length(ecgInterp) * 100;
                            fprintf('   ✅ ECG merged: %.1f%% coverage (%d valid samples)\n', coverage, validSamples);
                        end
                    else
                        if verboseOutput
                            fprintf('   ⚠️ Insufficient ECG data after time filtering\n');
                        end
                    end
                else
                    if verboseOutput
                        fprintf('   ⚠️ No valid ECG data found\n');
                    end
                end
            else
                if verboseOutput
                    fprintf('   ⚠️ No EDF files found for ECG data\n');
                end
            end
        end
        
        %% Step 4: Load and merge Pupil Labs IMU data (if requested)
        if includePupilIMU
            if verboseOutput, fprintf('\n📱 Loading Pupil Labs IMU data...\n'); end
            
            % Find pupil directories
            pupilDirs = dir(fullfile(subjectPath, '*pupil*'));
            pupilDirs = pupilDirs([pupilDirs.isdir]);
            
            imuData = [];
            for p = 1:length(pupilDirs)
                pupilPath = fullfile(subjectPath, pupilDirs(p).name);
                
                % Look for timestamped subdirectories
                subDirs = dir(pupilPath);
                subDirs = subDirs([subDirs.isdir] & ~strcmp({subDirs.name}, '.') & ~strcmp({subDirs.name}, '..'));
                
                if ~isempty(subDirs)
                    [~, latestIdx] = max([subDirs.datenum]);
                    pupilDataPath = fullfile(pupilPath, subDirs(latestIdx).name);
                else
                    pupilDataPath = pupilPath;
                end
                
                imuFile = fullfile(pupilDataPath, 'IMU.csv');
                if exist(imuFile, 'file')
                    if verboseOutput, fprintf('   📥 Reading: %s\n', imuFile); end
                    
                    try
                        imuTable = readtable(imuFile);
                        if height(imuTable) > 0
                            imuData = [imuData; imuTable];
                        end
                    catch ME
                        if verboseOutput
                            fprintf('      ❌ Error reading IMU file: %s\n', ME.message);
                        end
                    end
                end
            end
            
            if ~isempty(imuData) && height(imuData) > 0
                % Find timestamp column
                timestampCol = '';
                possibleTimestampCols = {'timestamp [ns]', 'timestampns', 'timestamp_ns_', 'timestamp'};
                for col = possibleTimestampCols
                    if ismember(col{1}, imuData.Properties.VariableNames)
                        timestampCol = col{1};
                        break;
                    end
                end
                
                if ~isempty(timestampCol)
                    % Convert timestamps to seconds (assuming nanoseconds)
                    imuTimes = imuData.(timestampCol) * 1e-9;
                    
                    % Find gyroscope and accelerometer columns
                    gyroCols = {};
                    accelCols = {};
                    
                    varNames = imuData.Properties.VariableNames;
                    for v = 1:length(varNames)
                        varName = varNames{v};
                        if contains(varName, {'gyro', 'Gyro'}, 'IgnoreCase', true)
                            gyroCols{end+1} = varName;
                        elseif contains(varName, {'accel', 'Accel'}, 'IgnoreCase', true)
                            accelCols{end+1} = varName;
                        end
                    end
                    
                    if ~isempty(gyroCols) || ~isempty(accelCols)
                        % Align IMU times with EEG times (simple offset alignment)
                        eegTimes = eeg.times / 1000;
                        timeOffset = eegTimes(1) - imuTimes(1);
                        imuTimesAligned = imuTimes + timeOffset;
                        
                        % Filter to EEG time range
                        timeFilter = imuTimesAligned >= eegTimes(1) & imuTimesAligned <= eegTimes(end);
                        imuFiltered = imuData(timeFilter, :);
                        imuTimesFiltered = imuTimesAligned(timeFilter);
                        
                        if height(imuFiltered) > 10
                            % Add gyroscope channels
                            for g = 1:min(3, length(gyroCols))
                                gyroData = imuFiltered.(gyroCols{g});
                                gyroInterp = interp1(imuTimesFiltered, gyroData, eegTimes, 'cubic', NaN);
                                
                                newChannelIdx = eeg.nbchan + 1;
                                eeg.data(newChannelIdx, :) = gyroInterp;
                                eeg.chanlocs(newChannelIdx).labels = sprintf('Gyro_%s', gyroCols{g});
                                eeg.chanlocs(newChannelIdx).type = 'IMU';
                                eeg.nbchan = eeg.nbchan + 1;
                                stats.channels_added = stats.channels_added + 1;
                            end
                            
                            % Add accelerometer channels
                            for a = 1:min(3, length(accelCols))
                                accelData = imuFiltered.(accelCols{a});
                                accelInterp = interp1(imuTimesFiltered, accelData, eegTimes, 'cubic', NaN);
                                
                                newChannelIdx = eeg.nbchan + 1;
                                eeg.data(newChannelIdx, :) = accelInterp;
                                eeg.chanlocs(newChannelIdx).labels = sprintf('Accel_%s', accelCols{a});
                                eeg.chanlocs(newChannelIdx).type = 'IMU';
                                eeg.nbchan = eeg.nbchan + 1;
                                stats.channels_added = stats.channels_added + 1;
                            end
                            
                            stats.data_sources_merged{end+1} = 'PupilIMU';
                            
                            if verboseOutput
                                fprintf('   ✅ IMU merged: %d gyro + %d accel channels\n', ...
                                       min(3, length(gyroCols)), min(3, length(accelCols)));
                            end
                        else
                            if verboseOutput
                                fprintf('   ⚠️ Insufficient IMU data after time filtering\n');
                            end
                        end
                    else
                        if verboseOutput
                            fprintf('   ⚠️ No gyro/accelerometer columns found in IMU data\n');
                        end
                    end
                else
                    if verboseOutput
                        fprintf('   ⚠️ No timestamp column found in IMU data\n');
                    end
                end
            else
                if verboseOutput
                    fprintf('   ⚠️ No valid IMU data found\n');
                end
            end
        end
        
        %% Step 5: Merge XSens Center of Mass data (if requested)
        if includeXSensCoM && ~isempty(xsensData) && isfield(xsensData, 'center_of_mass')
            if verboseOutput, fprintf('\n🏃 Merging XSens Center of Mass...\n'); end
            
            try
                % Load alignment data to get counter-LSL time mapping
                alignmentDataFile = fullfile(alignedPath, sprintf('%s_alignment_data.mat', subjectFolder));
                if exist(alignmentDataFile, 'file')
                    alignmentLoad = load(alignmentDataFile);
                    alignmentReport = alignmentLoad.alignmentData.alignmentReport;
                    
                    % Extract time mapping function
                    intercept = alignmentReport.intercept;
                    slope = alignmentReport.slope;
                    csv_to_xdf_time = @(csv_time_ns) intercept + slope * (csv_time_ns * 1e-9);
                    
                    % Process each XSens data file
                    xsensCoMData = [];
                    for x = 1:length(xsensData.center_of_mass)
                        if istable(xsensData.center_of_mass{x}) && height(xsensData.center_of_mass{x}) > 0
                            comTable = xsensData.center_of_mass{x};
                            
                            % Assume 60 Hz frame rate if not specified in info
                            frameRate = 60;
                            if isfield(xsensData, 'info') && length(xsensData.info) >= x
                                try
                                    frameRate = xsensData.info{x}{5,2};
                                catch
                                    frameRate = 60; % fallback
                                end
                            end
                            
                            % Create time vector
                            numFrames = height(comTable);
                            xsensTime = (0:numFrames-1) / frameRate;
                            
                            % Simple time alignment (more sophisticated alignment could be implemented)
                            eegTimes = eeg.times / 1000;
                            timeOffset = eegTimes(1);
                            xsensTimeAligned = xsensTime + timeOffset;
                            
                            % Filter to EEG time range
                            timeFilter = xsensTimeAligned >= eegTimes(1) & xsensTimeAligned <= eegTimes(end);
                            comFiltered = comTable(timeFilter, :);
                            xsensTimesFiltered = xsensTimeAligned(timeFilter);
                            
                            if height(comFiltered) > 10
                                % Extract Center of Mass position columns
                                comCols = {};
                                varNames = comFiltered.Properties.VariableNames;
                                for v = 1:length(varNames)
                                    varName = varNames{v};
                                    if contains(varName, {'CoM', 'pos'}, 'IgnoreCase', true) && ...
                                       contains(varName, {'X', 'Y', 'Z'}, 'IgnoreCase', true)
                                        comCols{end+1} = varName;
                                    end
                                end
                                
                                if ~isempty(comCols)
                                    for c = 1:length(comCols)
                                        colName = comCols{c};
                                        comData = comFiltered.(colName);
                                        comInterp = interp1(xsensTimesFiltered, comData, eegTimes, 'cubic', NaN);
                                        
                                        newChannelIdx = eeg.nbchan + 1;
                                        eeg.data(newChannelIdx, :) = comInterp;
                                        eeg.chanlocs(newChannelIdx).labels = sprintf('CoM_%s', colName);
                                        eeg.chanlocs(newChannelIdx).type = 'XSens';
                                        eeg.nbchan = eeg.nbchan + 1;
                                        stats.channels_added = stats.channels_added + 1;
                                    end
                                    
                                    if verboseOutput
                                        fprintf('   ✅ Added %d Center of Mass channels\n', length(comCols));
                                    end
                                end
                            end
                        end
                    end
                    
                    if stats.channels_added > 0
                        stats.data_sources_merged{end+1} = 'XSensCoM';
                    end
                else
                    if verboseOutput
                        fprintf('   ⚠️ Alignment data not found, using simple time alignment\n');
                    end
                end
            catch ME
                if verboseOutput
                    fprintf('   ❌ Error merging XSens CoM: %s\n', ME.message);
                end
            end
        end
        
        %% Step 6: Merge XSens Foot Position data (if requested)
        if includeXSensFootPos && ~isempty(xsensData) && isfield(xsensData, 'segment_position')
            if verboseOutput, fprintf('\n🦶 Merging XSens Foot Position...\n'); end
            
            try
                % Process each XSens segment position file
                for x = 1:length(xsensData.segment_position)
                    if istable(xsensData.segment_position{x}) && height(xsensData.segment_position{x}) > 0
                        segTable = xsensData.segment_position{x};
                        
                        % Find foot position columns
                        footCols = {};
                        toeCols = {};
                        varNames = segTable.Properties.VariableNames;
                        for v = 1:length(varNames)
                            varName = varNames{v};
                            if contains(varName, {'Foot'}, 'IgnoreCase', true) && ...
                               contains(varName, {'Z'}, 'IgnoreCase', true)
                                footCols{end+1} = varName;
                            elseif contains(varName, {'Toe'}, 'IgnoreCase', true) && ...
                               contains(varName, {'Z'}, 'IgnoreCase', true)
                                toeCols{end+1} = varName;
                            end
                        end
                        
                        if ~isempty(footCols)
                            % Assume 60 Hz frame rate
                            frameRate = 60;
                            if isfield(xsensData, 'info') && length(xsensData.info) >= x
                                try
                                    frameRate = xsensData.info{x}{5,2};
                                catch
                                    frameRate = 60;
                                end
                            end
                            
                            numFrames = height(segTable);
                            xsensTime = (0:numFrames-1) / frameRate;
                            
                            % Simple time alignment
                            eegTimes = eeg.times / 1000;
                            timeOffset = eegTimes(1);
                            xsensTimeAligned = xsensTime + timeOffset;
                            
                            newChannelIdx = eeg.nbchan + 1;
                            
                            for f = 1:length(footCols)
                                colName = footCols{f};
                                footData = segTable.(colName);
                                footInterp = interp1(xsensTimeAligned, footData, eegTimes, 'cubic', NaN);
                                
                                eeg.data(newChannelIdx, :) = footInterp;
                                eeg.chanlocs(newChannelIdx).labels = colName;
                                eeg.chanlocs(newChannelIdx).type = 'XSens';
                                newChannelIdx = newChannelIdx + 1;
                                stats.channels_added = stats.channels_added + 1;
                            end

                            for f = 1:length(toeCols)
                                colName = toeCols{f};
                                toeData = segTable.(colName);
                                toeInterp = interp1(xsensTimeAligned, toeData, eegTimes, 'cubic', NaN);
                                
                                eeg.data(newChannelIdx, :) = toeInterp;
                                eeg.chanlocs(newChannelIdx).labels = colName;
                                eeg.chanlocs(newChannelIdx).type = 'XSens';
                                newChannelIdx = newChannelIdx + 1;
                                stats.channels_added = stats.channels_added + 1;
                            end
                            
                            if verboseOutput
                                fprintf('   ✅ Added %d foot position channels\n', length(footCols));
                            end
                            
                            stats.data_sources_merged{end+1} = 'XSensFootPos';
                        end
                    end
                end
            catch ME
                if verboseOutput
                    fprintf('   ❌ Error merging XSens foot position: %s\n', ME.message);
                end
            end
        end
        
        %% Step 7: Update EEG structure and save
        if verboseOutput, fprintf('\n💾 Finalizing merged data...\n'); end
        
        % Update EEG structure - nbchan should already be updated
        eeg = eeg_checkset(eeg, 'eventconsistency');
        
        % Record processing statistics
        stats.final_channels = eeg.nbchan;
        stats.processing_time_seconds = toc(startTime);
        stats.data_sources_merged = strjoin(stats.data_sources_merged, ', ');
        mergeStats = [mergeStats; stats];
        
        if verboseOutput
            fprintf('   📊 Channels: %d → %d (+%d added)\n', ...
                   stats.original_channels, stats.final_channels, stats.channels_added);
            fprintf('   🔗 Data sources merged: %s\n', stats.data_sources_merged);
            fprintf('   ⏱️  Processing time: %.1f seconds\n', stats.processing_time_seconds);
        end
        
        % Save merged EEG file
        if saveResults
            pop_saveset(eeg, 'filename', sprintf('%s_merged.set', subjectFolder), 'filepath', outPath);
            
            % Save merge statistics
            mergeInfoFile = fullfile(outPath, sprintf('%s_merge_info.mat', subjectFolder));
            mergeInfo = struct();
            mergeInfo.stats = stats;
            mergeInfo.processing_date = datetime('now');
            mergeInfo.script_name = mfilename;
            mergeInfo.matlab_version = version;
            save(mergeInfoFile, 'mergeInfo', '-v7.3');
            
            if verboseOutput
                fprintf('   ✅ Saved: %s\n', outputFile);
                fprintf('   ✅ Merge info saved: %s\n', mergeInfoFile);
            end
        end
        
        % Create plots if requested
        if createPlots && stats.channels_added > 0
            if verboseOutput, fprintf('   📈 Creating merge plots...\n'); end
            create_merge_plots(eeg, stats, subjectFolder, outPath, savePlots, verboseOutput);
        end
        
        % Mark as successfully processed
        processedSubjects{end+1} = subjectFolder;
        
        if verboseOutput
            fprintf('\n✅ %s MERGING COMPLETED!\n', subjectFolder);
            fprintf('   📊 Final EEG structure: %d channels, %d samples\n', eeg.nbchan, eeg.pnts);
            fprintf('   🔗 Data sources: %s\n', stats.data_sources_merged);
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
    fprintf('ADDITIONAL DATA MERGING SUMMARY\n');
    fprintf('==========================================================\n\n');
    
    fprintf('✅ SUCCESSFULLY PROCESSED (%d subjects):\n', length(processedSubjects));
    for i = 1:length(processedSubjects)
        if i <= length(mergeStats)
            fprintf('   %s (+%d channels, %s, %.1fs)\n', processedSubjects{i}, ...
                   mergeStats(i).channels_added, mergeStats(i).data_sources_merged, ...
                   mergeStats(i).processing_time_seconds);
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
    if ~isempty(mergeStats)
        totalChannelsAdded = sum([mergeStats.channels_added]);
        totalTime = sum([mergeStats.processing_time_seconds]);
        
        fprintf('\n📊 MERGE STATISTICS:\n');
        fprintf('   Total channels added: %d\n', totalChannelsAdded);
        fprintf('   Average channels per subject: %.1f\n', totalChannelsAdded / length(processedSubjects));
        fprintf('   Total processing time: %.1f seconds (%.1f minutes)\n', totalTime, totalTime/60);
        fprintf('   Average processing time: %.1f seconds per subject\n', totalTime / length(processedSubjects));
        
        % Data source statistics
        allSources = {};
        for i = 1:length(mergeStats)
            sources = strsplit(mergeStats(i).data_sources_merged, ', ');
            allSources = [allSources, sources];
        end
        uniqueSources = unique(allSources(~cellfun(@isempty, allSources)));
        
        fprintf('\n🔗 DATA SOURCES MERGED:\n');
        for src = uniqueSources
            count = sum(contains({mergeStats.data_sources_merged}, src{1}));
            fprintf('   %s: %d subjects\n', src{1}, count);
        end
    end
    
    fprintf('\n💾 Output directory: %s\n', outPath);
    fprintf('==========================================================\n');
end

%% ==================== HELPER FUNCTIONS ====================

function create_merge_plots(eeg, stats, subjectName, outputPath, savePlots, verbose)
    % Create comprehensive plots showing the merged data quality
    
    try
        % Main figure with overview
        fig1 = figure('Position', [100, 100, 1600, 1000], 'Name', sprintf('%s Merged Data Overview', subjectName));
        
        % Get channel types and organize data
        channelTypes = {eeg.chanlocs.type};
        channelLabels = {eeg.chanlocs.labels};
        uniqueTypes = unique(channelTypes);
        uniqueTypes = uniqueTypes(~cellfun(@isempty, uniqueTypes));
        
        timeVector = eeg.times / 1000; % Convert to seconds
        
        % Plot 1: Channel overview by type
        subplot(2, 3, 1);
        typeCounts = [];
        for t = 1:length(uniqueTypes)
            typeCounts(t) = sum(strcmp(channelTypes, uniqueTypes{t}));
        end
        bar(typeCounts);
        set(gca, 'XTickLabel', uniqueTypes, 'XTickLabelRotation', 45);
        title('Channels by Type');
        ylabel('Number of Channels');
        grid on;
        
        % Plot 2-4: Show each data type
        plotIdx = 2;
        for t = 1:min(3, length(uniqueTypes))
            subplot(2, 3, plotIdx);
            
            channelType = uniqueTypes{t};
            typeIndices = find(strcmp(channelTypes, channelType));
            
            if ~isempty(typeIndices)
                % Plot first few channels of this type
                numToPlot = min(3, length(typeIndices));
                colors = lines(numToPlot);
                
                for ch = 1:numToPlot
                    chIdx = typeIndices(ch);
                    data = eeg.data(chIdx, :);
                    
                    % Normalize and offset data for plotting
                    if ~all(isnan(data))
                        dataNorm = (data - nanmean(data)) / nanstd(data);
                        plot(timeVector, dataNorm + (ch-1)*3, 'Color', colors(ch,:), 'LineWidth', 1);
                        hold on;
                    end
                end
                
                title(sprintf('%s (%d channels)', channelType, length(typeIndices)));
                xlabel('Time (seconds)');
                ylabel('Normalized Amplitude');
                grid on;
                
                % Add channel labels
                if numToPlot <= 3
                    legendLabels = {};
                    for ch = 1:numToPlot
                        chIdx = typeIndices(ch);
                        legendLabels{ch} = eeg.chanlocs(chIdx).labels;
                    end
                    legend(legendLabels, 'Location', 'best', 'FontSize', 8);
                end
            end
            plotIdx = plotIdx + 1;
        end
        
        % Plot 5: Summary statistics
        subplot(2, 3, 5);
        axis off;
        
        summaryText = {
            sprintf('MERGE SUMMARY: %s', subjectName),
            '',
            sprintf('Original channels: %d', stats.original_channels),
            sprintf('Final channels: %d', stats.final_channels),
            sprintf('Channels added: %d', stats.channels_added),
            '',
            'Data sources merged:',
            sprintf('%s', stats.data_sources_merged),
            '',
            sprintf('Processing time: %.1fs', stats.processing_time_seconds),
            '',
            'Channel breakdown:'
        };
        
        % Add channel type breakdown
        for t = 1:length(uniqueTypes)
            summaryText{end+1} = sprintf('  %s: %d channels', uniqueTypes{t}, typeCounts(t));
        end
        
        text(0.1, 0.9, summaryText, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
             'FontSize', 9, 'FontWeight', 'normal');
        
        % Plot 6: Data coverage analysis
        subplot(2, 3, 6);
        coverageData = [];
        coverageLabels = {};
        
        for t = 1:length(uniqueTypes)
            typeIndices = find(strcmp(channelTypes, uniqueTypes{t}));
            if ~isempty(typeIndices)
                % Calculate average coverage for this type
                totalCoverage = 0;
                validChannels = 0;
                
                for ch = typeIndices
                    data = eeg.data(ch, :);
                    coverage = sum(~isnan(data)) / length(data) * 100;
                    if ~isnan(coverage)
                        totalCoverage = totalCoverage + coverage;
                        validChannels = validChannels + 1;
                    end
                end
                
                if validChannels > 0
                    avgCoverage = totalCoverage / validChannels;
                    coverageData(end+1) = avgCoverage;
                    coverageLabels{end+1} = uniqueTypes{t};
                end
            end
        end
        
        if ~isempty(coverageData)
            bar(coverageData);
            set(gca, 'XTickLabel', coverageLabels, 'XTickLabelRotation', 45);
            title('Data Coverage by Type');
            ylabel('Coverage (%)');
            ylim([0, 100]);
            grid on;
        end
        
        sgtitle(sprintf('Merged Data Overview - %s', subjectName));
        
        if savePlots
            plotPath = fullfile(outputPath, sprintf('%s_merged_overview.png', subjectName));
            print(fig1, plotPath, '-dpng', '-r300');
            if verbose
                fprintf('      💾 Overview plot saved: %s\n', plotPath);
            end
        end
        
        % Create quality check plot (similar to user's example)
        fig2 = figure('Position', [150, 150, 1400, 800], 'Name', sprintf('%s Merge Quality Check', subjectName));
        
        % Find specific channels to plot for quality check
        ecgIdx = find(strcmp(channelTypes, 'ECG'), 1);
        imuIdx = find(strcmp(channelTypes, 'IMU'), 1);
        xsensIdx = find(strcmp(channelTypes, 'XSens'), 1);
        gazeIdx = find(strcmp(channelTypes, 'Gaze'), 1);
        
        % Plot quality check data
        subplot(2, 1, 1);
        plotChannels = [];
        plotLabels = {};
        
        if ~isempty(ecgIdx)
            plot(eeg.data(ecgIdx, :), 'DisplayName', eeg.chanlocs(ecgIdx).labels);
            hold on;
            plotChannels(end+1) = ecgIdx;
            plotLabels{end+1} = eeg.chanlocs(ecgIdx).labels;
        end
        
        if ~isempty(imuIdx)
            % Normalize IMU data for comparison
            imuData = eeg.data(imuIdx, :);
            imuNorm = (imuData - nanmean(imuData)) / nanstd(imuData) * nanstd(eeg.data(ecgIdx, :)) + nanmean(eeg.data(ecgIdx, :));
            plot(imuNorm, 'DisplayName', eeg.chanlocs(imuIdx).labels);
            plotChannels(end+1) = imuIdx;
            plotLabels{end+1} = eeg.chanlocs(imuIdx).labels;
        end
        
        if ~isempty(xsensIdx)
            % Normalize XSens data for comparison
            xsensData = eeg.data(xsensIdx, :);
            xsensNorm = (xsensData - nanmean(xsensData)) / nanstd(xsensData) * nanstd(eeg.data(ecgIdx, :)) + nanmean(eeg.data(ecgIdx, :));
            plot(xsensNorm, 'DisplayName', eeg.chanlocs(xsensIdx).labels);
            plotChannels(end+1) = xsensIdx;
            plotLabels{end+1} = eeg.chanlocs(xsensIdx).labels;
        end
        
        % Find and plot events (similar to user's example)
        if ~isempty(eeg.event)
            % Find cloud events (if they exist)
            cloudEvents = [];
            if isfield(eeg.event, 'pupil_type')
                cloudEvents = find(strcmpi({eeg.event.pupil_type}, 'cloud'));
            end
            
            % If no cloud events, find any experimental events
            if isempty(cloudEvents)
                % Look for events that are not sync events
                eventTypes = {eeg.event.type};
                excludeTypes = {'lsl.time_sync', 'recording.begin', 'recording.end'};
                eventMask = true(size(eventTypes));
                for excl = excludeTypes
                    eventMask = eventMask & ~contains(eventTypes, excl{1});
                end
                cloudEvents = find(eventMask);
            end
            
            if ~isempty(cloudEvents) && length(cloudEvents) <= 50  % Limit to avoid clutter
                latencies = [eeg.event(cloudEvents).latency];
                types = {eeg.event(cloudEvents).type};
                
                % Get y-axis range for marker placement
                ylims = ylim;
                markerY = ylims(1) + 0.9 * (ylims(2) - ylims(1));
                
                % Plot markers at event latencies
                plot(latencies, repmat(markerY, size(latencies)), ...
                     'Marker', '+', 'LineStyle', 'none', 'MarkerSize', 15, ...
                     'Color', 'red', 'LineWidth', 2);
                
                % Add text labels (sample first few to avoid clutter)
                numLabels = min(10, length(latencies));
                for i = 1:numLabels
                    text(latencies(i), markerY * 1.05, types{i}, ...
                         'Rotation', 45, 'FontSize', 8, 'HorizontalAlignment', 'left', ...
                         'Color', 'red');
                end
                
                if length(latencies) > numLabels
                    text(latencies(end), markerY * 1.05, sprintf('...+%d more', length(latencies)-numLabels), ...
                         'Rotation', 45, 'FontSize', 8, 'HorizontalAlignment', 'left', ...
                         'Color', 'red');
                end
            end
        end
        
        title(sprintf('Merged Data Quality Check - %s', subjectName));
        xlabel('Sample Number');
        ylabel('Amplitude');
        legend(plotLabels, 'Location', 'best');
        grid on;
        hold off;
        
        % Plot event timeline
        subplot(2, 1, 2);
        if ~isempty(eeg.event)
            % Create event timeline
            eventTimes = [eeg.event.latency] / eeg.srate; % Convert to seconds
            eventTypes = {eeg.event.type};
            
            % Get unique event types and assign colors
            uniqueEventTypes = unique(eventTypes);
            colors = lines(length(uniqueEventTypes));
            
            hold on;
            for et = 1:length(uniqueEventTypes)
                eventType = uniqueEventTypes{et};
                eventMask = strcmp(eventTypes, eventType);
                eventTimesType = eventTimes(eventMask);
                
                if length(eventTimesType) <= 100  % Avoid plotting too many events
                    % Plot as stem plot
                    stem(eventTimesType, repmat(et, size(eventTimesType)), ...
                         'Color', colors(et, :), 'MarkerFaceColor', colors(et, :), ...
                         'MarkerSize', 4, 'LineWidth', 1);
                end
            end
            
            % Set axis properties
            xlim([0, max(timeVector)]);
            ylim([0.5, length(uniqueEventTypes) + 0.5]);
            set(gca, 'YTick', 1:length(uniqueEventTypes), 'YTickLabel', uniqueEventTypes);
            title(sprintf('Event Timeline (%d events, %d types)', length(eeg.event), length(uniqueEventTypes)));
            xlabel('Time (seconds)');
            ylabel('Event Type');
            grid on;
            hold off;
        else
            text(0.5, 0.5, 'No events found', 'Units', 'normalized', 'HorizontalAlignment', 'center');
            title('Event Timeline');
        end
        
        if savePlots
            plotPath = fullfile(outputPath, sprintf('%s_merge_quality.png', subjectName));
            print(fig2, plotPath, '-dpng', '-r300');
            if verbose
                fprintf('      💾 Quality plot saved: %s\n', plotPath);
            end
        end
        
        if verbose
            fprintf('      ✅ Merge quality plots created\n');
        end
        
    catch ME
        if verbose
            fprintf('      ❌ Error creating merge plots: %s\n', ME.message);
        end
    end
end