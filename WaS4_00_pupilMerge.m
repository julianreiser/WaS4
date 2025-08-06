%% EEG-PUPIL LABS PREPROCESSING PIPELINE - ENHANCED
% One-step alignment and merging script for EEG and Pupil Labs data
% Run this as the first step in your preprocessing pipeline
%
% This script:
% 1. Finds all subjects in the data directory
% 2. Loads EEG (BDF), XDF, and Pupil Labs CSV data (with multi-file support)
% 3. Performs time alignment using LSL sync events
% 4. Merges all data into aligned EEG files
% 5. Saves results and generates reports
%
% ENHANCED FEATURES:
% - Multiple XDF file support with timestamp reset detection
% - Multiple Pupil Labs folder support
% - Robust error handling and fallback options
% - Better file merging logic

clear; clc; close all;

dataPath = '/Volumes/ergo/rohdaten/RohDaten_GRAIL/WaS4/DATA/';  % Main data directory
outputPath = '/Volumes/Work4TB/Seafile/WaS4/data/aligned/';  % Output directory for aligned files

%% EEG-PUPIL LABS PREPROCESSING PIPELINE

% Processing options
processAllSubjects = false;  % Set to false to process specific subjects
specificSubjects = [46];  % Only used if processAllSubjects = false
overwriteExisting = false;  % Set to true to reprocess existing files

% Quality thresholds
minSyncEvents = 2;  % Minimum number of sync events required (reduced to 2 for recording.begin/end)
maxRMSE = 0.1;  % Maximum acceptable RMSE in seconds (much stricter for sync)
minRSquared = 0.99;  % Minimum acceptable R-squared (should be very high for good sync)

% Output options
saveResults = true;
createPlots = true;
savePlots = true;
verboseOutput = true;

%% ==================== INITIALIZATION ====================

if verboseOutput
    fprintf('==========================================================\n');
    fprintf('EEG-PUPIL LABS PREPROCESSING PIPELINE (ENHANCED)\n');
    fprintf('==========================================================\n\n');
    fprintf('📁 Data path: %s\n', dataPath);
    fprintf('💾 Output path: %s\n\n', outputPath);
end

% Create output directory
if ~exist(outputPath, 'dir')
    mkdir(outputPath);
    if verboseOutput, fprintf('📁 Created output directory: %s\n\n', outputPath); end
end

% Find all subject folders
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
alignmentQuality = [];

%% ==================== MAIN PROCESSING LOOP ====================

for s = 1:length(subjectList)
    subjectFolder = subjectList{s};
    subjectPath = fullfile(dataPath, subjectFolder);
    
    if verboseOutput
        fprintf('==========================================================\n');
        fprintf('🔄 PROCESSING: %s (%d/%d)\n', subjectFolder, s, length(subjectList));
        fprintf('==========================================================\n\n');
    end
    
    try
        % Check if already processed
        outputFile = fullfile(outputPath, sprintf('%s_aligned.set', subjectFolder));
        if exist(outputFile, 'file') && ~overwriteExisting
            if verboseOutput
                fprintf('⏭️  Already processed (use overwriteExisting=true to reprocess)\n\n');
            end
            skippedSubjects{end+1} = subjectFolder;
            continue;
        end
        
        %% Step 1: Load EEG data (ENHANCED - Multiple BDF support)
        if verboseOutput, fprintf('🧠 Loading EEG data...\n'); end
        
        eegFiles = dir(fullfile(subjectPath, '*.bdf'));
        if isempty(eegFiles)
            error('No BDF files found');
        end
        
        % Load and merge EEG files
        eegData = {};
        ampCNTOffsets = [0];  % Track counter offsets for multi-file support
        
        for f = 1:length(eegFiles)
            filepath = fullfile(subjectPath, eegFiles(f).name);
            if verboseOutput, fprintf('   📥 Reading %s\n', eegFiles(f).name); end
            eegData{f} = pop_biosig(filepath);
        end
        
        % Merge EEG files if multiple
        eeg = eegData{1};
        cntIdx = find(strcmp({eeg.chanlocs.labels}, 'CNT'));
        if isempty(cntIdx)
            error('No CNT (counter) channel found in EEG data');
        end
        
        if length(eegData) > 1
            if verboseOutput, fprintf('   🔗 Merging %d EEG files...\n', length(eegData)); end
            for f = 2:length(eegData)
                newTime = eeg.times(end) + eegData{f}.times + (1000/eeg.srate);
                eeg.times = [eeg.times newTime];
                
                % Adjust counter value to continue from the last file
                eegData{f}.data(cntIdx, :) = eegData{f}.data(cntIdx, :) + eeg.data(cntIdx,end);
                ampCNTOffsets = [ampCNTOffsets eeg.data(cntIdx, end)]; % Store offsets for LSL adjustment
                
                eeg.data = [eeg.data, eegData{f}.data];
                eeg.pnts = eeg.pnts + eegData{f}.pnts;
            end
        end
        
        if verboseOutput
            fprintf('   ✅ EEG loaded: %d channels, %d samples, %.1f Hz (%d files merged)\n\n', ...
                   eeg.nbchan, eeg.pnts, eeg.srate, length(eegFiles));
        end
        
        %% Step 2: Load XDF data (ENHANCED - Multiple XDF support with timestamp reset detection)
        if verboseOutput, fprintf('📊 Loading XDF data...\n'); end
        
        xdfFiles = dir(fullfile(subjectPath, '*.xdf'));
        if isempty(xdfFiles)
            error('No XDF files found');
        end
        
        % Load all XDF files
        streamsTmp = {};
        for f = 1:length(xdfFiles)
            filepath = fullfile(subjectPath, xdfFiles(f).name);
            if verboseOutput, fprintf('   📥 Reading %s\n', xdfFiles(f).name); end
            streamsTmp{f} = load_xdf(filepath);
        end
        
        % Handle multiple XDF files with timestamp reset detection
        if length(streamsTmp) > 1
            if verboseOutput, fprintf('   🔧 Detecting timestamp resets and merging %d XDF files...\n', length(streamsTmp)); end
            
            for f = 2:length(streamsTmp)
                reset = false;
                
                % Find counter/sync streams in consecutive files
                idx1 = cellfun(@(x) contains(x.info.name, {'Counter', 'Sync', 'PROX'}, 'IgnoreCase', true), streamsTmp{f-1});
                idx2 = cellfun(@(x) contains(x.info.name, {'Counter', 'Sync', 'PROX'}, 'IgnoreCase', true), streamsTmp{f});
                
                if any(idx1) && any(idx2)
                    if streamsTmp{f}{idx2}.time_stamps(1) < streamsTmp{f-1}{idx1}.time_stamps(end)
                        reset = true;
                        if verboseOutput, fprintf('   ⚠️ Timestamp reset detected between files %d and %d\n', f-1, f); end
                    end
                end
                
                % Fix timestamps if reset detected
                if reset
                    if verboseOutput, fprintf('   🔧 Correcting timestamp reset...\n'); end
                    
                    % Calculate time offset using EEG counter if possible
                    if length(eegData) > 1
                        % Use EEG counter-based correction
                        idxLast = find(idx1, 1);
                        idxCurrent = find(idx2, 1);
                        
                        if ~isempty(idxLast) && ~isempty(idxCurrent)
                            % Find corresponding times in EEG data
                            lastCounterVal = streamsTmp{f-1}{idxLast}.time_series(end);
                            currentCounterVal = streamsTmp{f}{idxCurrent}.time_series(1);
                            
                            % Find these counter values in EEG
                            lastEEGIdx = find(eeg.data(cntIdx,:) == lastCounterVal, 1, 'last');
                            currentEEGIdx = find(eeg.data(cntIdx,:) == currentCounterVal + ampCNTOffsets(f), 1, 'first');
                            
                            if ~isempty(lastEEGIdx) && ~isempty(currentEEGIdx)
                                timeDiffEEG = (eeg.times(currentEEGIdx) - eeg.times(lastEEGIdx)) / 1000; % Convert to seconds
                                lslTimeOffset = streamsTmp{f-1}{idxLast}.time_stamps(end) + timeDiffEEG - streamsTmp{f}{idxCurrent}.time_stamps(1);
                            else
                                % Fallback: use simple continuation
                                lslTimeOffset = streamsTmp{f-1}{idxLast}.time_stamps(end) - streamsTmp{f}{idxCurrent}.time_stamps(1) + 0.001;
                                if verboseOutput, fprintf('   ⚠️ Using fallback timestamp correction\n'); end
                            end
                        else
                            error('Cannot find counter streams for timestamp correction');
                        end
                    else
                        % Simple fallback for single EEG file
                        lslTimeOffset = streamsTmp{f-1}{idx1}.time_stamps(end) - streamsTmp{f}{idx2}.time_stamps(1) + 0.001;
                        if verboseOutput, fprintf('   ⚠️ Using simple timestamp correction\n'); end
                    end
                    
                    % Apply offset to all streams in current and subsequent files
                    for nf = f:length(streamsTmp)
                        for d = 1:length(streamsTmp{nf})
                            if ~isempty(streamsTmp{nf}{d}.time_stamps)
                                streamsTmp{nf}{d}.time_stamps = streamsTmp{nf}{d}.time_stamps + lslTimeOffset;
                            end
                        end
                    end
                    
                    if verboseOutput, fprintf('   ✅ Applied time offset: %.6f seconds\n', lslTimeOffset); end
                end
            end
        end
        
        % Merge XDF streams from all files
        streams = streamsTmp{1};
        if length(streamsTmp) > 1
            if verboseOutput, fprintf('   🔗 Merging XDF streams...\n'); end
            
            for f = 2:length(streamsTmp)
                nextStream = streamsTmp{f};
                for n = 1:length(nextStream)
                    % Find matching stream in main streams
                    matchFound = false;
                    for id = 1:length(streams)
                        if strcmp(streams{id}.info.name, nextStream{n}.info.name)
                            % Handle counter adjustment for EEG streams
                            if strcmp(streams{id}.info.type, "EEG") && ~isempty(nextStream{n}.time_series)
                                if nextStream{n}.time_series(1) < streams{id}.time_series(end)
                                    % Adjust EEG counter values
                                    offsetIdx = find(ampCNTOffsets + nextStream{n}.time_series(1) > streams{id}.time_series(end), 1, 'first');
                                    if ~isempty(offsetIdx)
                                        nextStream{n}.time_series = nextStream{n}.time_series + ampCNTOffsets(offsetIdx);
                                    end
                                end
                            end
                            
                            % Merge time stamps and data
                            streams{id}.time_stamps = [streams{id}.time_stamps nextStream{n}.time_stamps];
                            streams{id}.time_series = [streams{id}.time_series nextStream{n}.time_series];
                            matchFound = true;
                            break;
                        end
                    end
                    
                    % If no matching stream found, add as new stream
                    if ~matchFound
                        streams{end+1} = nextStream{n};
                        if verboseOutput, fprintf('     + Added new stream: %s\n', nextStream{n}.info.name); end
                    end
                end
            end
        end
        
        if verboseOutput
            fprintf('   ✅ XDF loaded: %d streams (%d files merged)\n\n', length(streams), length(xdfFiles));
        end
        
        %% Step 3: Load Pupil Labs CSV data (ENHANCED - Multiple folder support)
        if verboseOutput, fprintf('👁️ Loading Pupil Labs data...\n'); end
        
        % Find all pupil cloud folders (enhanced pattern matching)
        pupilDirs = dir(fullfile(subjectPath, '*pupil*'));
        pupilDirs = pupilDirs([pupilDirs.isdir]);
        
        if isempty(pupilDirs)
            error('No Pupil Labs folders found');
        end
        
        if verboseOutput, fprintf('   🔍 Found %d Pupil Labs folders\n', length(pupilDirs)); end
        
        % Process all pupil directories
        pupilData = struct();
        pupilData.events = table();
        pupilData.gaze = table();
        pupilData.fixations = table();
        pupilData.blinks = table();
        
        for pupilIdx = 1:length(pupilDirs)
            pupilDir = fullfile(subjectPath, pupilDirs(pupilIdx).name);
            if verboseOutput, fprintf('   📂 Processing folder: %s\n', pupilDirs(pupilIdx).name); end
            
            % Find the actual data folder (timestamped subfolder or direct folder)
            dataFolders = dir(pupilDir);
            dataFolders = dataFolders([dataFolders.isdir] & ~strcmp({dataFolders.name}, '.') & ~strcmp({dataFolders.name}, '..'));
            
            if ~isempty(dataFolders)
                [~, latest_idx] = max([dataFolders.datenum]);
                pupilDataPath = fullfile(pupilDir, dataFolders(latest_idx).name);
                if verboseOutput, fprintf('     📁 Using subfolder: %s\n', dataFolders(latest_idx).name); end
            else
                pupilDataPath = pupilDir;
                if verboseOutput, fprintf('     📁 Using direct folder\n'); end
            end
            
            % Load CSV files with enhanced file detection
            csvFiles = {'events.csv', 'gaze.csv', 'time_aligned_gaze.csv', 'fixations.csv', 'blinks.csv'};
            currentPupilData = struct();
            
            for i = 1:length(csvFiles)
                csvPath = fullfile(pupilDataPath, csvFiles{i});
                if exist(csvPath, 'file')
                    try
                        % Check file size first
                        fileInfo = dir(csvPath);
                        if fileInfo.bytes < 100
                            if verboseOutput
                                fprintf('     ⚠️ %s too small (%d bytes), skipping\n', csvFiles{i}, fileInfo.bytes);
                            end
                            continue;
                        end
                        
                        data = readtable(csvPath);
                        if height(data) > 0
                            fieldName = csvFiles{i}(1:end-4);
                            % Handle special case for time_aligned_gaze
                            if strcmp(fieldName, 'time_aligned_gaze')
                                fieldName = 'gaze';
                            end
                            currentPupilData.(fieldName) = data;
                            if verboseOutput
                                fprintf('     📥 Loaded %s: %d rows\n', csvFiles{i}, height(data));
                            end
                        else
                            if verboseOutput
                                fprintf('     ⚠️ %s is empty\n', csvFiles{i});
                            end
                        end
                    catch ME
                        if verboseOutput
                            fprintf('     ❌ Could not load %s: %s\n', csvFiles{i}, ME.message);
                        end
                    end
                else
                    if verboseOutput && i <= 2 % Only warn for critical files
                        fprintf('     ⚠️ File not found: %s\n', csvFiles{i});
                    end
                end
            end
            
            % Merge data from current folder into main pupilData structure
            csvDataTypes = {'events', 'gaze', 'fixations', 'blinks'};
            for dataType = csvDataTypes
                dataTypeStr = dataType{1};
                if isfield(currentPupilData, dataTypeStr)
                    if height(pupilData.(dataTypeStr)) == 0
                        pupilData.(dataTypeStr) = currentPupilData.(dataTypeStr);
                    else
                        % Concatenate with existing data
                        try
                            pupilData.(dataTypeStr) = vertcat(pupilData.(dataTypeStr), currentPupilData.(dataTypeStr));
                            if verboseOutput
                                fprintf('     🔗 Merged %s data (%d new rows)\n', dataTypeStr, height(currentPupilData.(dataTypeStr)));
                            end
                        catch ME
                            if verboseOutput
                                fprintf('     ⚠️ Could not merge %s data: %s\n', dataTypeStr, ME.message);
                            end
                        end
                    end
                end
            end
        end
        
        % Check if we have required data
        if ~isfield(pupilData, 'events') || height(pupilData.events) == 0
            error('events.csv is required for time synchronization but not found or empty');
        end
        
        % Summary of loaded Pupil Labs data
        if verboseOutput
            fprintf('   📊 PUPIL LABS DATA SUMMARY:\n');
            dataTypes = {'events', 'gaze', 'fixations', 'blinks'};
            for dataType = dataTypes
                if isfield(pupilData, dataType{1}) && height(pupilData.(dataType{1})) > 0
                    fprintf('     - %s: %d rows\n', dataType{1}, height(pupilData.(dataType{1})));
                else
                    fprintf('     - %s: not available\n', dataType{1});
                end
            end
            fprintf('\n');
        end
        
        %% Step 4: Find LSL time sync events and perform alignment
        if verboseOutput, fprintf('🔄 Performing time alignment...\n'); end
        
        % Extract sync events from XDF (try multiple methods)
        xdfSyncEvents = [];
        xdfSyncTimes = [];
        syncMethod = '';
        
        for i = 1:length(streams)
            if isfield(streams{i}, 'time_series') && isfield(streams{i}, 'info') && ...
               (isfield(streams{i}.info, 'type') && ...
                (strcmpi(streams{i}.info.type, 'Markers') || strcmpi(streams{i}.info.type, 'Event')))
                
                events = streams{i}.time_series;
                times = streams{i}.time_stamps;
                
                if verboseOutput
                    fprintf('   🔍 Checking stream "%s" (type: %s)\n', ...
                           streams{i}.info.name, streams{i}.info.type);
                end
                
                if iscell(events)
                    % Method 1: Look for lsl.time_sync events (preferred)
                    syncMask = contains(events, 'lsl.time_sync');
                    if any(syncMask)
                        xdfSyncEvents = [xdfSyncEvents; events(syncMask)];
                        xdfSyncTimes = [xdfSyncTimes; times(syncMask)'];
                        if verboseOutput
                            fprintf('      ✅ Found %d lsl.time_sync events\n', sum(syncMask));
                        end
                    end
                    
                    % Method 2: Look for recording events (fallback)
                    recordingMask = contains(events, 'recording.begin') | contains(events, 'recording.end');
                    if any(recordingMask)
                        tempEvents = events(recordingMask);
                        tempTimes = times(recordingMask)';
                        if verboseOutput
                            fprintf('      ✅ Found %d recording.begin/end events\n', sum(recordingMask));
                        end
                        
                        % If we don't have enough lsl.time_sync events, use recording events
                        if length(xdfSyncEvents) < minSyncEvents
                            xdfSyncEvents = [xdfSyncEvents; tempEvents];
                            xdfSyncTimes = [xdfSyncTimes; tempTimes];
                        end
                    end
                end
            end
        end
        
        % Determine which sync method we're using
        lslSyncCount = sum(contains(xdfSyncEvents, 'lsl.time_sync'));
        recordingCount = sum(contains(xdfSyncEvents, 'recording.'));
        
        if lslSyncCount >= minSyncEvents
            syncMethod = 'lsl.time_sync';
            if verboseOutput, fprintf('   🎯 Using LSL time sync events (%d found)\n', lslSyncCount); end
        elseif recordingCount >= 2
            syncMethod = 'recording.begin/end';
            if verboseOutput, fprintf('   🎯 Using recording begin/end events (%d found)\n', recordingCount); end
        elseif (lslSyncCount + recordingCount) >= minSyncEvents
            syncMethod = 'mixed';
            if verboseOutput, fprintf('   🎯 Using mixed sync method (%d LSL + %d recording)\n', lslSyncCount, recordingCount); end
        else
            syncMethod = 'insufficient';
        end
        
        % Extract corresponding sync events from CSV
        csvEvents = pupilData.events;
        csvSyncEvents = [];
        csvSyncTimesNs = [];
        
        % Look for the same types of events in CSV as we found in XDF
        if contains(syncMethod, 'lsl.time_sync') || strcmp(syncMethod, 'mixed')
            csvSyncMask = contains(csvEvents.name, 'lsl.time_sync');
            if any(csvSyncMask)
                csvSyncEvents = [csvSyncEvents; csvEvents.name(csvSyncMask)];
                
                % Find timestamp column
                timestampCol = find_timestamp_column(csvEvents);
                csvSyncTimesNs = [csvSyncTimesNs; csvEvents.(timestampCol)(csvSyncMask)];
            end
        end
        
        if contains(syncMethod, 'recording') || strcmp(syncMethod, 'mixed') || strcmp(syncMethod, 'insufficient')
            csvRecordingMask = contains(csvEvents.name, 'recording.begin') | contains(csvEvents.name, 'recording.end');
            if any(csvRecordingMask)
                csvSyncEvents = [csvSyncEvents; csvEvents.name(csvRecordingMask)];
                
                if isempty(csvSyncTimesNs)
                    timestampCol = find_timestamp_column(csvEvents);
                end
                csvSyncTimesNs = [csvSyncTimesNs; csvEvents.(timestampCol)(csvRecordingMask)];
            end
        end
        
        if verboseOutput
            fprintf('   📊 Found %d sync events in XDF\n', length(xdfSyncEvents));
            fprintf('   📊 Found %d sync events in CSV\n', length(csvSyncEvents));
            
            if ~isempty(xdfSyncTimes) && ~isempty(csvSyncTimesNs)
                fprintf('   📊 XDF time range: %.3f - %.3f seconds\n', min(xdfSyncTimes), max(xdfSyncTimes));
                fprintf('   📊 CSV time range: %.0f - %.0f ns (%.3f - %.3f s since epoch)\n', ...
                       min(csvSyncTimesNs), max(csvSyncTimesNs), ...
                       min(csvSyncTimesNs)*1e-9, max(csvSyncTimesNs)*1e-9);
            end
        end
        
        % Check minimum requirements
        if length(xdfSyncEvents) < 2 || length(csvSyncEvents) < 2
            error('Insufficient sync events (need at least 2, found XDF:%d, CSV:%d). Method: %s', ...
                  length(xdfSyncEvents), length(csvSyncEvents), syncMethod);
        end
        
        % Sort events by time
        [xdfSyncTimes, xdfSortIdx] = sort(xdfSyncTimes);
        xdfSyncEvents = xdfSyncEvents(xdfSortIdx);
        
        [csvSyncTimesNs, csvSortIdx] = sort(csvSyncTimesNs);
        csvSyncEvents = csvSyncEvents(csvSortIdx);
        
        % Enhanced event matching based on sync method
        if strcmp(syncMethod, 'lsl.time_sync')
            [matchedXdf, matchedCsvNs] = match_sync_events_robust(xdfSyncEvents, xdfSyncTimes, ...
                                                                 csvSyncEvents, csvSyncTimesNs, verboseOutput);
        elseif strcmp(syncMethod, 'recording.begin/end')
            [matchedXdf, matchedCsvNs] = match_recording_events(xdfSyncEvents, xdfSyncTimes, ...
                                                              csvSyncEvents, csvSyncTimesNs, verboseOutput);
        else % mixed method
            [matchedXdf, matchedCsvNs] = match_mixed_events(xdfSyncEvents, xdfSyncTimes, ...
                                                          csvSyncEvents, csvSyncTimesNs, verboseOutput);
        end
        
        if length(matchedXdf) < 2
            error('Could not match sufficient sync events (need at least 2, found %d)', length(matchedXdf));
        end
        
        % Convert CSV times to seconds for the linear fit
        matchedCsvSeconds = matchedCsvNs * 1e-9;
        
        % PROPER TIME ALIGNMENT METHOD
        % The key insight: we need to work with time differences, not absolute times
        % because CSV uses Unix timestamps and XDF uses LSL timestamps
        
        if length(matchedXdf) >= 2
            % Method 1: Use time differences between consecutive sync events
            xdfDiffs = diff(matchedXdf);
            csvDiffs = diff(matchedCsvSeconds);
            
            if verboseOutput
                fprintf('   🔧 Using time differences method for alignment\n');
                fprintf('      XDF time diffs: %.3f ± %.3f seconds\n', mean(xdfDiffs), std(xdfDiffs));
                fprintf('      CSV time diffs: %.3f ± %.3f seconds\n', mean(csvDiffs), std(csvDiffs));
            end
            
            % Calculate slope (should be close to 1.0 for good sync)
            if length(xdfDiffs) > 1
                slope = xdfDiffs \ csvDiffs;  % Robust slope estimation
            else
                slope = xdfDiffs / csvDiffs;
            end
            
            % Calculate intercept using the first matched pair
            intercept = matchedXdf(1) - slope * matchedCsvSeconds(1);
            
            % Quality metrics using differences
            if length(xdfDiffs) > 1
                predicted_diffs = slope * csvDiffs;
                residuals_diff = xdfDiffs - predicted_diffs;
                rmse_diff = sqrt(mean(residuals_diff.^2));
                r_squared_diff = 1 - sum(residuals_diff.^2) / sum((xdfDiffs - mean(xdfDiffs)).^2);
                
                % Overall quality using predicted vs actual times
                predicted_xdf = intercept + slope * matchedCsvSeconds;
                residuals = matchedXdf - predicted_xdf;
                rmse = sqrt(mean(residuals.^2));
                r_squared = 1 - sum(residuals.^2) / sum((matchedXdf - mean(matchedXdf)).^2);
                
                if verboseOutput
                    fprintf('      Time differences - RMSE: %.6f s, R²: %.6f\n', rmse_diff, r_squared_diff);
                end
            else
                % Fallback for only 2 sync events
                predicted_xdf = intercept + slope * matchedCsvSeconds;
                residuals = matchedXdf - predicted_xdf;
                rmse = sqrt(mean(residuals.^2));
                r_squared = 1 - sum(residuals.^2) / sum((matchedXdf - mean(matchedXdf)).^2);
            end
            
        else
            error('Need at least 2 matched sync events for time alignment');
        end
        
        % Validate the alignment makes sense
        if abs(slope - 1.0) > 0.1
            warning('Clock drift detected: slope = %.6f (should be close to 1.0)', slope);
        end
        
        % Adjust quality thresholds based on sync method
        if strcmp(syncMethod, 'lsl.time_sync')
            rmseThreshold = maxRMSE;
            r2Threshold = 0.99;
        else
            % Be more lenient with recording.begin/end events
            rmseThreshold = maxRMSE * 5;  % Allow up to 0.5s RMSE
            r2Threshold = 0.95;           % Lower R² threshold
        end
        
        if rmse > rmseThreshold
            warning('High RMSE (%.3f > %.3f) - alignment may be poor (method: %s)', rmse, rmseThreshold, syncMethod);
        end
        
        if r_squared < r2Threshold
            warning('Low R-squared (%.6f < %.6f) - alignment may be poor (method: %s)', r_squared, r2Threshold, syncMethod);
        end
        
        alignmentReport = struct();
        alignmentReport.subject = subjectFolder;
        alignmentReport.intercept = intercept;
        alignmentReport.slope = slope;
        alignmentReport.r_squared = r_squared;
        alignmentReport.rmse = rmse;
        alignmentReport.sync_events_used = length(matchedXdf);
        alignmentReport.matched_pairs_xdf = matchedXdf;
        alignmentReport.matched_pairs_csv_ns = matchedCsvNs;
        alignmentReport.matched_pairs_csv_s = matchedCsvSeconds;
        alignmentReport.timestamp_column = timestampCol;
        alignmentReport.clock_drift = abs(slope - 1.0);
        alignmentReport.sync_method = syncMethod;
        alignmentReport.files_processed.eeg_files = length(eegFiles);
        alignmentReport.files_processed.xdf_files = length(xdfFiles);
        alignmentReport.files_processed.pupil_folders = length(pupilDirs);
        
        alignmentQuality = [alignmentQuality; alignmentReport];
        
        if verboseOutput
            fprintf('   ✅ Time mapping (%s): intercept=%.6f, slope=%.9f\n', syncMethod, intercept, slope);
            fprintf('   📊 Quality: R²=%.6f, RMSE=%.6f sec (%d events)\n', r_squared, rmse, length(matchedXdf));
            fprintf('   ⏱️  Clock drift: %.2f%% (slope deviation from 1.0)\n', abs(slope-1.0)*100);
            fprintf('   📁 Files: %d EEG, %d XDF, %d Pupil folders\n', length(eegFiles), length(xdfFiles), length(pupilDirs));
            fprintf('\n');
        end
        
        %% Step 5: Apply alignment to all Pupil Labs data
        if verboseOutput, fprintf('🔧 Applying time alignment to Pupil Labs data...\n'); end
        
        % Function to convert CSV timestamps (nanoseconds) to XDF time (seconds)
        csv_to_xdf_time = @(csv_time_ns) intercept + slope * (csv_time_ns * 1e-9);
        
        % Align events
        pupilData.events.timestamp_xdf = csv_to_xdf_time(pupilData.events.(timestampCol));
        if verboseOutput, fprintf('   ✅ Aligned %d events\n', height(pupilData.events)); end
        
        % Align gaze data
        if isfield(pupilData, 'gaze') && height(pupilData.gaze) > 0
            gazeTimestampCol = find_timestamp_column(pupilData.gaze);
            if ~isempty(gazeTimestampCol)
                pupilData.gaze.timestamp_xdf = csv_to_xdf_time(pupilData.gaze.(gazeTimestampCol));
                if verboseOutput, fprintf('   ✅ Aligned %d gaze samples\n', height(pupilData.gaze)); end
            else
                if verboseOutput, fprintf('   ⚠️ No timestamp column found in gaze data\n'); end
            end
        end
        
        % Align fixations
        if isfield(pupilData, 'fixations') && height(pupilData.fixations) > 0
            start_col = find(contains(pupilData.fixations.Properties.VariableNames, 'startTimestamp', 'IgnoreCase', true));
            end_col = find(contains(pupilData.fixations.Properties.VariableNames, 'endTimestamp', 'IgnoreCase', true));
            
            if ~isempty(start_col)
                col_name = pupilData.fixations.Properties.VariableNames{start_col(1)};
                pupilData.fixations.startTimestamp_xdf = csv_to_xdf_time(pupilData.fixations.(col_name));
            end
            if ~isempty(end_col)
                col_name = pupilData.fixations.Properties.VariableNames{end_col(1)};
                pupilData.fixations.endTimestamp_xdf = csv_to_xdf_time(pupilData.fixations.(col_name));
            end
            if verboseOutput, fprintf('   ✅ Aligned %d fixations\n', height(pupilData.fixations)); end
        end
        
        % Align blinks
        if isfield(pupilData, 'blinks') && height(pupilData.blinks) > 0
            start_col = find(contains(pupilData.blinks.Properties.VariableNames, 'startTimestamp', 'IgnoreCase', true));
            end_col = find(contains(pupilData.blinks.Properties.VariableNames, 'endTimestamp', 'IgnoreCase', true));
            
            if ~isempty(start_col)
                col_name = pupilData.blinks.Properties.VariableNames{start_col(1)};
                pupilData.blinks.startTimestamp_xdf = csv_to_xdf_time(pupilData.blinks.(col_name));
            end
            if ~isempty(end_col)
                col_name = pupilData.blinks.Properties.VariableNames{end_col(1)};
                pupilData.blinks.endTimestamp_xdf = csv_to_xdf_time(pupilData.blinks.(col_name));
            end
            if verboseOutput, fprintf('   ✅ Aligned %d blinks\n', height(pupilData.blinks)); end
        end
        
        if verboseOutput, fprintf('\n'); end
        
        %% Step 6: Merge aligned Pupil Labs data into EEG structure
        if verboseOutput, fprintf('🔀 Merging data into EEG structure...\n'); end
        
        % Find EEG counter stream in XDF
        eegIdx = [];
        for i = 1:length(streams)
            if isfield(streams{i}.info, 'name') && ...
               (contains(streams{i}.info.name, 'Counter', 'IgnoreCase', true) || ...
                contains(streams{i}.info.name, 'Sync', 'IgnoreCase', true) || ...
                contains(streams{i}.info.name, 'PROX', 'IgnoreCase', true))
                eegIdx = i;
                break;
            end
        end
        
        if isempty(eegIdx)
            error('Could not find EEG counter stream in XDF data');
        end
        
        % Match counter values between EEG and XDF
        [~, idxAmp, idxLSL] = intersect(eeg.data(cntIdx,:), streams{eegIdx}.time_series);
        counterXdfTimes = streams{eegIdx}.time_stamps(idxLSL);
        
        % Define EEG recording time bounds (used for gaze filtering and plotting)
        eegTimeStart = min(counterXdfTimes);
        eegTimeEnd = max(counterXdfTimes);
        
        if verboseOutput
            fprintf('   📊 Found %d matching counter values\n', length(counterXdfTimes));
            fprintf('   📊 EEG time range: %.3f - %.3f seconds (%.1f min)\n', ...
                   eegTimeStart, eegTimeEnd, (eegTimeEnd - eegTimeStart)/60);
        end
        
        % Add gaze data to EEG
        if isfield(pupilData, 'gaze') && height(pupilData.gaze) > 0 && ismember('timestamp_xdf', pupilData.gaze.Properties.VariableNames)
            gaze_x_col = find_column_by_names(pupilData.gaze.Properties.VariableNames, {'gazeX', 'gaze_x', 'gazeX_px', 'gazeXpx'});
            gaze_y_col = find_column_by_names(pupilData.gaze.Properties.VariableNames, {'gazeY', 'gaze_y', 'gazeY_px', 'gazeYpx'});
            
            if ~isempty(gaze_x_col) && ~isempty(gaze_y_col)
                if verboseOutput
                    fprintf('   📊 Gaze time range: %.3f - %.3f seconds\n', ...
                           min(pupilData.gaze.timestamp_xdf), max(pupilData.gaze.timestamp_xdf));
                end
                
                % Filter gaze data to EEG time range with small buffer
                timeBuffer = 1.0; % 1 second buffer on each side
                gazeTimeFilter = (pupilData.gaze.timestamp_xdf >= (eegTimeStart - timeBuffer)) & ...
                                (pupilData.gaze.timestamp_xdf <= (eegTimeEnd + timeBuffer));
                
                % Also remove invalid gaze points
                validGaze = gazeTimeFilter & ...
                           ~isnan(pupilData.gaze.(gaze_x_col)) & ...
                           ~isnan(pupilData.gaze.(gaze_y_col)) & ...
                           ~isnan(pupilData.gaze.timestamp_xdf);
                
                if verboseOutput
                    fprintf('   📊 Gaze samples: %d total, %d in time range, %d valid\n', ...
                           height(pupilData.gaze), sum(gazeTimeFilter), sum(validGaze));
                end
                
                if sum(validGaze) > 10  % Need at least some valid data
                    % Extract filtered gaze data
                    gazeTimesValid = pupilData.gaze.timestamp_xdf(validGaze);
                    gazeXValid = pupilData.gaze.(gaze_x_col)(validGaze);
                    gazeYValid = pupilData.gaze.(gaze_y_col)(validGaze);
                    
                    % Check for overlap between gaze and EEG times
                    if max(gazeTimesValid) < eegTimeStart || min(gazeTimesValid) > eegTimeEnd
                        if verboseOutput
                            fprintf('   ❌ No temporal overlap between gaze and EEG data\n');
                        end
                    else
                        % Interpolate gaze data to EEG counter times
                        gazeX = interp1(gazeTimesValid, gazeXValid, counterXdfTimes, 'linear', NaN);
                        gazeY = interp1(gazeTimesValid, gazeYValid, counterXdfTimes, 'linear', NaN);
                        
                        % Check interpolation quality
                        validInterp = ~isnan(gazeX) & ~isnan(gazeY);
                        interpCoverage = sum(validInterp) / length(validInterp) * 100;
                        
                        if verboseOutput
                            fprintf('   📊 Interpolation coverage: %.1f%% (%d/%d samples)\n', ...
                                   interpCoverage, sum(validInterp), length(validInterp));
                        end
                        
                        if interpCoverage > 1.0  % At least 1% coverage
                            % Add channels to EEG
                            newIdx = eeg.nbchan + 1;
                            eeg.data(newIdx:newIdx+1, :) = NaN(2, size(eeg.data, 2));
                            eeg.data(newIdx:newIdx+1, idxAmp) = [gazeX; gazeY];
                            
                            eeg.chanlocs(newIdx).labels = 'GazeX';
                            eeg.chanlocs(newIdx).type = 'Gaze';
                            eeg.chanlocs(newIdx+1).labels = 'GazeY';
                            eeg.chanlocs(newIdx+1).type = 'Gaze';
                            
                            eeg.nbchan = eeg.nbchan + 2;
                            
                            if verboseOutput
                                fprintf('   ✅ Added gaze X/Y channels to EEG (%.1f%% coverage)\n', interpCoverage);
                            end
                        else
                            if verboseOutput
                                fprintf('   ⚠️ Poor interpolation coverage (%.1f%%), skipping gaze data\n', interpCoverage);
                            end
                        end
                    end
                else
                    if verboseOutput
                        fprintf('   ⚠️ Not enough valid gaze data in EEG time range (%d samples)\n', sum(validGaze));
                    end
                end
            else
                if verboseOutput
                    fprintf('   ⚠️ Could not find gaze X/Y columns. Available: %s\n', ...
                           strjoin(pupilData.gaze.Properties.VariableNames, ', '));
                end
            end
        else
            if verboseOutput, fprintf('   ⚠️ No gaze data or timestamp_xdf column found\n'); end
        end
        
        % Add events to EEG
        eventsAdded = 0;
        blinksAdded = 0;
        fixationsAdded = 0;
        saccadesAdded = 0;
        
        if isfield(pupilData, 'events') && ismember('timestamp_xdf', pupilData.events.Properties.VariableNames)
            % Convert XDF times to EEG sample indices
            eegTimes = eeg.times / 1000; % Convert to seconds
            eegTimes = eegTimes - eegTimes(idxAmp(1)) + counterXdfTimes(1); % Align with XDF time
            
            eventCount = length(eeg.event);
            eventsAdded = 0;
            
            % Add events from events.csv
            for i = 1:height(pupilData.events)
                [timeDiff, sampleIdx] = min(abs(eegTimes - pupilData.events.timestamp_xdf(i)));
                
                % Only add events that are within reasonable time bounds
                if timeDiff < 0.1  % Within 100ms
                    eventCount = eventCount + 1;
                    eeg.event(eventCount).latency = sampleIdx;
                    eeg.event(eventCount).duration = 1;
                    eeg.event(eventCount).type = char(pupilData.events.name(i));
                    eeg.event(eventCount).source = 'pupil_events';
                    
                    if ismember('type', pupilData.events.Properties.VariableNames)
                        eeg.event(eventCount).pupil_type = char(pupilData.events.type(i));
                    end
                    
                    eventsAdded = eventsAdded + 1;
                end
            end
            
            if verboseOutput, fprintf('   ✅ Added %d events from events.csv to EEG\n', eventsAdded); end
        end
        
        % Add blinks as events
        if isfield(pupilData, 'blinks') && height(pupilData.blinks) > 0 && ismember('startTimestamp_xdf', pupilData.blinks.Properties.VariableNames)
            if verboseOutput, fprintf('   🔄 Adding blink events to EEG...\n'); end
            
            blinksAdded = 0;
            for i = 1:height(pupilData.blinks)
                % Use start timestamp for event timing
                [timeDiff, sampleIdx] = min(abs(eegTimes - pupilData.blinks.startTimestamp_xdf(i)));
                
                if timeDiff < 0.1  % Within 100ms
                    eventCount = eventCount + 1;
                    eeg.event(eventCount).latency = sampleIdx;
                    eeg.event(eventCount).type = 'blink';
                    eeg.event(eventCount).source = 'pupil_blinks';
                    
                    % Add blink-specific information
                    if ismember('duration_ms_', pupilData.blinks.Properties.VariableNames)
                        duration_ms = pupilData.blinks.duration_ms_(i);
                        if ~isnan(duration_ms)
                            eeg.event(eventCount).duration = round(duration_ms * eeg.srate / 1000); % Convert to samples
                            eeg.event(eventCount).blink_duration_ms = duration_ms;
                        else
                            eeg.event(eventCount).duration = 1;
                        end
                    elseif ismember('durationms', pupilData.blinks.Properties.VariableNames)
                        duration_ms = pupilData.blinks.durationms(i);
                        if ~isnan(duration_ms)
                            eeg.event(eventCount).duration = round(duration_ms * eeg.srate / 1000);
                            eeg.event(eventCount).blink_duration_ms = duration_ms;
                        else
                            eeg.event(eventCount).duration = 1;
                        end
                    else
                        eeg.event(eventCount).duration = 1;
                    end
                    
                    % Add start and end timestamps if available
                    if ismember('endTimestamp_xdf', pupilData.blinks.Properties.VariableNames)
                        eeg.event(eventCount).blink_end_time_xdf = pupilData.blinks.endTimestamp_xdf(i);
                    end
                    eeg.event(eventCount).blink_start_time_xdf = pupilData.blinks.startTimestamp_xdf(i);
                    
                    blinksAdded = blinksAdded + 1;
                end
            end
            
            if verboseOutput, fprintf('   ✅ Added %d blink events to EEG\n', blinksAdded); end
        end
        
        % Add fixations as events
        if isfield(pupilData, 'fixations') && height(pupilData.fixations) > 0 && ismember('startTimestamp_xdf', pupilData.fixations.Properties.VariableNames)
            if verboseOutput, fprintf('   🔄 Adding fixation events to EEG...\n'); end
            
            fixationsAdded = 0;
            for i = 1:height(pupilData.fixations)
                % Use start timestamp for event timing
                [timeDiff, sampleIdx] = min(abs(eegTimes - pupilData.fixations.startTimestamp_xdf(i)));
                
                if timeDiff < 0.1  % Within 100ms
                    eventCount = eventCount + 1;
                    eeg.event(eventCount).latency = sampleIdx;
                    eeg.event(eventCount).type = 'fixation';
                    eeg.event(eventCount).source = 'pupil_fixations';
                    
                    % Add fixation-specific information with flexible column name matching
                    duration_cols = {'duration_ms_', 'durationms', 'duration'};
                    for col = duration_cols
                        if ismember(col{1}, pupilData.fixations.Properties.VariableNames)
                            duration_ms = pupilData.fixations.(col{1})(i);
                            if ~isnan(duration_ms)
                                eeg.event(eventCount).duration = round(duration_ms * eeg.srate / 1000);
                                eeg.event(eventCount).fixation_duration_ms = duration_ms;
                            else
                                eeg.event(eventCount).duration = 1;
                            end
                            break;
                        end
                    end
                    if ~isfield(eeg.event(eventCount), 'duration')
                        eeg.event(eventCount).duration = 1;
                    end
                    
                    % Add fixation coordinates with flexible naming
                    x_cols = {'fixation_x_px_', 'fixationX_px_', 'fixationXpx', 'fixation_x'};
                    y_cols = {'fixation_y_px_', 'fixationY_px_', 'fixationYpx', 'fixation_y'};
                    
                    for col = x_cols
                        if ismember(col{1}, pupilData.fixations.Properties.VariableNames)
                            eeg.event(eventCount).fixation_x_px = pupilData.fixations.(col{1})(i);
                            break;
                        end
                    end
                    
                    for col = y_cols
                        if ismember(col{1}, pupilData.fixations.Properties.VariableNames)
                            eeg.event(eventCount).fixation_y_px = pupilData.fixations.(col{1})(i);
                            break;
                        end
                    end
                    
                    % Add gaze angles if available
                    if ismember('azimuth_deg_', pupilData.fixations.Properties.VariableNames)
                        eeg.event(eventCount).fixation_azimuth_deg = pupilData.fixations.azimuth_deg_(i);
                    elseif ismember('azimuthdeg', pupilData.fixations.Properties.VariableNames)
                        eeg.event(eventCount).fixation_azimuth_deg = pupilData.fixations.azimuthdeg(i);
                    end
                    
                    if ismember('elevation_deg_', pupilData.fixations.Properties.VariableNames)
                        eeg.event(eventCount).fixation_elevation_deg = pupilData.fixations.elevation_deg_(i);
                    elseif ismember('elevationdeg', pupilData.fixations.Properties.VariableNames)
                        eeg.event(eventCount).fixation_elevation_deg = pupilData.fixations.elevationdeg(i);
                    end
                    
                    % Add fixation ID if available
                    id_cols = {'fixation_id', 'fixationId'};
                    for col = id_cols
                        if ismember(col{1}, pupilData.fixations.Properties.VariableNames)
                            eeg.event(eventCount).fixation_id = pupilData.fixations.(col{1})(i);
                            break;
                        end
                    end
                    
                    % Add timestamps
                    eeg.event(eventCount).fixation_start_time_xdf = pupilData.fixations.startTimestamp_xdf(i);
                    if ismember('endTimestamp_xdf', pupilData.fixations.Properties.VariableNames)
                        eeg.event(eventCount).fixation_end_time_xdf = pupilData.fixations.endTimestamp_xdf(i);
                    end
                    
                    fixationsAdded = fixationsAdded + 1;
                end
            end
            
            if verboseOutput, fprintf('   ✅ Added %d fixation events to EEG\n', fixationsAdded); end
        end
        
        % Add saccades as events (if available)
        if isfield(pupilData, 'saccades') && height(pupilData.saccades) > 0 && ismember('startTimestamp_xdf', pupilData.saccades.Properties.VariableNames)
            if verboseOutput, fprintf('   🔄 Adding saccade events to EEG...\n'); end
            
            for i = 1:height(pupilData.saccades)
                % Use start timestamp for event timing
                [timeDiff, sampleIdx] = min(abs(eegTimes - pupilData.saccades.startTimestamp_xdf(i)));
                
                if timeDiff < 0.1  % Within 100ms
                    eventCount = eventCount + 1;
                    eeg.event(eventCount).latency = sampleIdx;
                    eeg.event(eventCount).type = 'saccade';
                    eeg.event(eventCount).source = 'pupil_saccades';
                    
                    % Add saccade-specific information
                    if ismember('duration_ms_', pupilData.saccades.Properties.VariableNames)
                        duration_ms = pupilData.saccades.duration_ms_(i);
                        if ~isnan(duration_ms)
                            eeg.event(eventCount).duration = round(duration_ms * eeg.srate / 1000); % Convert to samples
                            eeg.event(eventCount).saccade_duration_ms = duration_ms;
                        else
                            eeg.event(eventCount).duration = 1;
                        end
                    else
                        eeg.event(eventCount).duration = 1;
                    end
                    
                    % Add saccade amplitude and velocity
                    if ismember('amplitude_px_', pupilData.saccades.Properties.VariableNames)
                        eeg.event(eventCount).saccade_amplitude_px = pupilData.saccades.amplitude_px_(i);
                    end
                    
                    if ismember('peak_velocity_px_s_', pupilData.saccades.Properties.VariableNames)
                        eeg.event(eventCount).saccade_peak_velocity = pupilData.saccades.peak_velocity_px_s_(i);
                    end
                    
                    if ismember('mean_velocity_px_s_', pupilData.saccades.Properties.VariableNames)
                        eeg.event(eventCount).saccade_mean_velocity = pupilData.saccades.mean_velocity_px_s_(i);
                    end
                    
                    % Add timestamps
                    eeg.event(eventCount).saccade_start_time_xdf = pupilData.saccades.startTimestamp_xdf(i);
                    if ismember('endTimestamp_xdf', pupilData.saccades.Properties.VariableNames)
                        eeg.event(eventCount).saccade_end_time_xdf = pupilData.saccades.endTimestamp_xdf(i);
                    end
                    
                    saccadesAdded = saccadesAdded + 1;
                end
            end
            
            if verboseOutput, fprintf('   ✅ Added %d saccade events to EEG\n', saccadesAdded); end
        else
            if verboseOutput, fprintf('   ⚠️ No saccades data available\n'); end
        end
        
        % Final summary of all events added
        totalEventsAdded = eventsAdded + blinksAdded + fixationsAdded + saccadesAdded;
        if verboseOutput
            fprintf('   📊 TOTAL EVENTS ADDED: %d (Events:%d, Blinks:%d, Fixations:%d, Saccades:%d)\n', ...
                   totalEventsAdded, eventsAdded, blinksAdded, fixationsAdded, saccadesAdded);
        end
        
        % Update EEG structure
        eeg = eeg_checkset(eeg);
        if verboseOutput, fprintf('   ✅ EEG structure updated: %d channels total\n\n', eeg.nbchan); end
        
        %% Step 7: Create alignment plots
        if createPlots
            if verboseOutput, fprintf('📈 Creating alignment plots...\n'); end
            
            % Extract and analyze special events FIRST for debugging
            if verboseOutput, fprintf('   🔍 Analyzing special events...\n'); end
            
            % Show all events for debugging
            if verboseOutput
                fprintf('   📊 Total events in CSV: %d\n', height(pupilData.events));
                uniqueEventNames = unique(pupilData.events.name);
                fprintf('   📊 Unique event types: %d\n', length(uniqueEventNames));
                fprintf('   📋 All event types found:\n');
                for i = 1:length(uniqueEventNames)
                    count = sum(strcmp(pupilData.events.name, uniqueEventNames{i}));
                    fprintf('      - %s (%d)\n', uniqueEventNames{i}, count);
                end
            end
            
            % Extract special events with debugging
            specialEvents = extract_special_events(pupilData.events, verboseOutput);
            
            % CREATE STANDALONE SPECIAL EVENTS FIGURE FIRST
            if ~isempty(specialEvents) && height(specialEvents) > 0
                if verboseOutput, fprintf('   🎨 Creating special events figure...\n'); end
                
                try
                    % Convert special events to XDF time
                    specialEventsXdf = specialEvents;
                    specialEventsXdf.timestamp_xdf = csv_to_xdf_time(specialEvents.(timestampCol));
                    
                    % Create standalone special events figure
                    figSpecial = figure('Position', [200, 200, 1400, 600], 'Name', sprintf('%s - Special Events', subjectFolder));
                    
                    % Filter events to EEG time range
                    timeFilter = specialEventsXdf.timestamp_xdf >= eegTimeStart & specialEventsXdf.timestamp_xdf <= eegTimeEnd;
                    eventsInRange = specialEventsXdf(timeFilter, :);
                    
                    if height(eventsInRange) > 0
                        % Option 1: Make timeline plot larger and add text summary
                        subplot(1, 1, 1);  % Use full figure for timeline
                        create_simple_events_timeline(eventsInRange, eegTimeStart, eegTimeEnd, verboseOutput);
                        title(sprintf('Special Events Timeline (%d events)', height(eventsInRange)));
                        xlabel('Time (seconds since LSL start)');
                        ylabel('Events');
                        
                        % Add text summary on the plot
                        [eventCounts, eventTypes] = count_events_by_type(eventsInRange);
                        
                        % Create a neat text summary
                        eventSummary = {};
                        for i = 1:length(eventTypes)
                            eventSummary{end+1} = sprintf('%s (%d)', eventTypes{i}, eventCounts(i));
                        end
                        
                        % Add summary text box
                        textStr = sprintf('Events Found:\n%s', strjoin(eventSummary, '\n'));
                        
                        % Position text box in upper right
                        text(0.98, 0.98, textStr, 'Units', 'normalized', ...
                             'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
                             'BackgroundColor', 'white', 'EdgeColor', 'black', ...
                             'FontSize', 9, 'Interpreter', 'none');
                        
                        if savePlots
                            specialPlotPath = fullfile(outputPath, sprintf('%s_special_events.png', subjectFolder));
                            print(figSpecial, specialPlotPath, '-dpng', '-r300');
                            if verboseOutput, fprintf('   💾 Special events plot saved: %s\n', specialPlotPath); end
                        end
                        
                        if verboseOutput, fprintf('   ✅ Special events figure created with %d events\n', height(eventsInRange)); end
                    else
                        text(0.5, 0.5, sprintf('No special events in EEG time range\n(Found %d events outside range)', height(specialEventsXdf)), ...
                             'Units', 'normalized', 'HorizontalAlignment', 'center');
                        title('Special Events Timeline');
                        if verboseOutput, fprintf('   ⚠️ No special events in EEG time range\n'); end
                    end
                    
                catch ME
                    if verboseOutput, fprintf('   ❌ Error creating special events plot: %s\n', ME.message); end
                end
            else
                if verboseOutput, fprintf('   ⚠️ No special events found to plot\n'); end
            end
            
            % Create main alignment figure
            fig = figure('Position', [100, 100, 1400, 1000], 'Name', sprintf('%s Alignment Results', subjectFolder));
            
            % Plot 1: Sync events alignment
            subplot(2, 3, 1);
            scatter(alignmentReport.matched_pairs_csv_s, alignmentReport.matched_pairs_xdf, 100, 'ro', 'filled');
            hold on;
            csv_range = [min(alignmentReport.matched_pairs_csv_s), max(alignmentReport.matched_pairs_csv_s)];
            xdf_fitted = intercept + slope * csv_range;
            plot(csv_range, xdf_fitted, 'b-', 'LineWidth', 2);
            xlabel('CSV Time (seconds since Unix epoch)');
            ylabel('XDF Time (seconds since LSL start)');
            title(sprintf('Time Synchronization (R² = %.6f)', r_squared));
            grid on;
            legend('Sync Events', 'Fitted Line', 'Location', 'best');
            
            % Plot 2: Residuals
            subplot(2, 3, 2);
            plot(alignmentReport.matched_pairs_csv_s, residuals, 'ro-');
            xlabel('CSV Time (seconds since Unix epoch)');
            ylabel('Residual (seconds)');
            title(sprintf('Alignment Residuals (RMSE = %.6f s)', rmse));
            grid on;
            
            % Plot 3: Data overview
            subplot(2, 3, 3);
            dataTypes = {};
            dataCounts = [];
            
            fields = fieldnames(pupilData);
            for i = 1:length(fields)
                if istable(pupilData.(fields{i}))
                    dataTypes{end+1} = fields{i};
                    dataCounts(end+1) = height(pupilData.(fields{i}));
                end
            end
            
            if ~isempty(dataCounts)
                bar(dataCounts);
                set(gca, 'XTickLabel', dataTypes);
                xlabel('Data Type');
                ylabel('Number of Records');
                title('Pupil Labs Data Overview');
                grid on;
            end
            
            % Plot 4: Special Events Summary
            subplot(2, 3, 4);
            if ~isempty(specialEvents) && height(specialEvents) > 0
                [eventCounts, eventTypes] = count_events_by_type(specialEvents);
                
                if length(eventTypes) <= 8
                    bar(eventCounts);
                    set(gca, 'XTickLabel', eventTypes, 'XTickLabelRotation', 45);
                    title(sprintf('Special Events (%d total)', sum(eventCounts)));
                    ylabel('Count');
                else
                    % Show top events only
                    topN = min(8, length(eventTypes));
                    bar(eventCounts(1:topN));
                    set(gca, 'XTickLabel', eventTypes(1:topN), 'XTickLabelRotation', 45);
                    title(sprintf('Top Special Events (%d/%d types)', topN, length(eventTypes)));
                    ylabel('Count');
                end
            else
                text(0.5, 0.5, 'No special events found', 'Units', 'normalized', 'HorizontalAlignment', 'center');
                title('Special Events');
            end
            
            % Plot 5: Event Timeline (simple version)
            subplot(2, 3, 5);
            if ~isempty(specialEvents) && height(specialEvents) > 0
                try
                    specialEventsXdf = specialEvents;
                    specialEventsXdf.timestamp_xdf = csv_to_xdf_time(specialEvents.(timestampCol));
                    
                    % Simple scatter plot of events over time
                    timeFilter = specialEventsXdf.timestamp_xdf >= eegTimeStart & specialEventsXdf.timestamp_xdf <= eegTimeEnd;
                    eventsInRange = specialEventsXdf(timeFilter, :);
                    
                    if height(eventsInRange) > 0
                        scatter(eventsInRange.timestamp_xdf, ones(height(eventsInRange), 1), 50, 'filled');
                        xlim([eegTimeStart, eegTimeEnd]);
                        ylim([0.5, 1.5]);
                        xlabel('Time (seconds since LSL start)');
                        title(sprintf('Events Timeline (%d events)', height(eventsInRange)));
                        set(gca, 'YTick', []);
                        grid on;
                    else
                        text(0.5, 0.5, 'No events in EEG range', 'Units', 'normalized', 'HorizontalAlignment', 'center');
                        title('Events Timeline');
                    end
                catch
                    text(0.5, 0.5, 'Error plotting timeline', 'Units', 'normalized', 'HorizontalAlignment', 'center');
                    title('Events Timeline');
                end
            else
                text(0.5, 0.5, 'No events to plot', 'Units', 'normalized', 'HorizontalAlignment', 'center');
                title('Events Timeline');
            end
            
            % Plot 6: Summary text
            subplot(2, 3, 6);
            axis off;
            
            summaryText = {
                sprintf('ALIGNMENT SUMMARY: %s', subjectFolder),
                '',
                'TIME MAPPING:',
                sprintf('  Method: %s', alignmentReport.sync_method),
                sprintf('  Slope: %.9f', slope),
                sprintf('  Intercept: %.6f s', intercept),
                '',
                'QUALITY METRICS:',
                sprintf('  R-squared: %.6f', r_squared),
                sprintf('  RMSE: %.6f seconds', rmse),
                sprintf('  Sync events used: %d', length(matchedXdf)),
                '',
                'FILES PROCESSED:',
                sprintf('  EEG files: %d', alignmentReport.files_processed.eeg_files),
                sprintf('  XDF files: %d', alignmentReport.files_processed.xdf_files),
                sprintf('  Pupil folders: %d', alignmentReport.files_processed.pupil_folders),
                '',
                'DATA MERGED:',
                sprintf('  EEG channels: %d', eeg.nbchan),
                sprintf('  Total events in EEG: %d', length(eeg.event)),
                sprintf('    - Pupil events: %d', eventsAdded),
                sprintf('    - Blinks: %d', blinksAdded),
                sprintf('    - Fixations: %d', fixationsAdded),
                sprintf('    - Saccades: %d', saccadesAdded),
                '',
                'SPECIAL EVENTS:'
            };
            
            % Add special event summary
            if ~isempty(specialEvents) && height(specialEvents) > 0
                summaryText{end+1} = sprintf('  Total special events: %d', height(specialEvents));
                [eventCounts, eventTypes] = count_events_by_type(specialEvents);
                for i = 1:min(5, length(eventTypes))  % Show top 5 event types
                    summaryText{end+1} = sprintf('  %s: %d', eventTypes{i}, eventCounts(i));
                end
                if length(eventTypes) > 5
                    summaryText{end+1} = sprintf('  ... and %d more types', length(eventTypes) - 5);
                end
            else
                summaryText{end+1} = '  No special events found';
            end
            
            % Color code based on quality
            if r_squared > 0.99 && rmse < 0.1
                textColor = [0, 0.6, 0];  % Green
            elseif r_squared > 0.95 && rmse < 0.5
                textColor = [0.8, 0.4, 0];  % Orange
            else
                textColor = [0.8, 0, 0];  % Red
            end
            
            text(0.05, 0.95, summaryText, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
                 'FontSize', 9, 'FontWeight', 'normal', 'Color', textColor);
            
            if savePlots
                plotPath = fullfile(outputPath, sprintf('%s_alignment_plots.png', subjectFolder));
                print(fig, plotPath, '-dpng', '-r300');
                if verboseOutput, fprintf('   💾 Main plots saved: %s\n', plotPath); end
            end
            
            % CREATE ADDITIONAL DETAILED SPECIAL EVENTS FIGURE
            if ~isempty(specialEvents) && height(specialEvents) > 0
                if verboseOutput, fprintf('   🎨 Creating detailed special events figure...\n'); end
                create_detailed_special_events_figure(specialEvents, csv_to_xdf_time, ...
                                                    eegTimeStart, eegTimeEnd, subjectFolder, ...
                                                    outputPath, savePlots, verboseOutput);
            end
            
            if verboseOutput, fprintf('   ✅ Alignment plots created\n\n'); end
        end
        
        % re-order events
        eeg = eeg_checkset(eeg,'eventconsistency');
        
        %% Step 8: Save results
        if saveResults
            if verboseOutput, fprintf('💾 Saving results...\n'); end
            
            % Save aligned EEG file
            outputFile = fullfile(outputPath, sprintf('%s_aligned.set', subjectFolder));
            pop_saveset(eeg, 'filename', sprintf('%s_aligned.set', subjectFolder), 'filepath', outputPath);
            
            % Save alignment report and pupil data
            reportPath = fullfile(outputPath, sprintf('%s_alignment_data.mat', subjectFolder));
            alignmentData = struct();
            alignmentData.alignmentReport = alignmentReport;
            alignmentData.pupilData = pupilData;
            alignmentData.originalStreams = streams;
            alignmentData.processingInfo.timestamp = datetime('now');
            alignmentData.processingInfo.matlabVersion = version;
            alignmentData.processingInfo.script = mfilename;
            save(reportPath, 'alignmentData', '-v7.3');
            
            if verboseOutput
                fprintf('   ✅ EEG data saved: %s\n', outputFile);
                fprintf('   ✅ Alignment data saved: %s\n', reportPath);
            end
        end
        
        % Mark as successfully processed
        processedSubjects{end+1} = subjectFolder;
        
        if verboseOutput
            fprintf('\n✅ %s COMPLETED SUCCESSFULLY!\n', subjectFolder);
            fprintf('   Quality (R²): %.6f, RMSE: %.3f sec\n', r_squared, rmse);
            fprintf('   Files processed: %d EEG, %d XDF, %d Pupil folders\n', ...
                   alignmentReport.files_processed.eeg_files, ...
                   alignmentReport.files_processed.xdf_files, ...
                   alignmentReport.files_processed.pupil_folders);
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
        
        % close all plots
        close all
    end
end

%% ==================== FINAL SUMMARY ====================

if verboseOutput
    fprintf('\n==========================================================\n');
    fprintf('PREPROCESSING PIPELINE SUMMARY (ENHANCED)\n');
    fprintf('==========================================================\n\n');
    
    fprintf('✅ SUCCESSFULLY PROCESSED (%d subjects):\n', length(processedSubjects));
    for i = 1:length(processedSubjects)
        if i <= length(alignmentQuality)
            fprintf('   %s (R²=%.4f, RMSE=%.3fs, Files: %dE/%dX/%dP)\n', processedSubjects{i}, ...
                   alignmentQuality(i).r_squared, alignmentQuality(i).rmse, ...
                   alignmentQuality(i).files_processed.eeg_files, ...
                   alignmentQuality(i).files_processed.xdf_files, ...
                   alignmentQuality(i).files_processed.pupil_folders);
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
    if ~isempty(alignmentQuality)
        allRSquared = [alignmentQuality.r_squared];
        allRMSE = [alignmentQuality.rmse];
        
        fprintf('\n📊 QUALITY STATISTICS:\n');
        fprintf('   R² - Mean: %.4f, Std: %.4f, Range: %.4f-%.4f\n', ...
               mean(allRSquared), std(allRSquared), min(allRSquared), max(allRSquared));
        fprintf('   RMSE - Mean: %.3fs, Std: %.3fs, Range: %.3f-%.3fs\n', ...
               mean(allRMSE), std(allRMSE), min(allRMSE), max(allRMSE));
               
        % Multi-file statistics
        allEEGFiles = [alignmentQuality.files_processed];
        eegCounts = [allEEGFiles.eeg_files];
        xdfCounts = [allEEGFiles.xdf_files];
        pupilCounts = [allEEGFiles.pupil_folders];
        
        fprintf('\n📁 MULTI-FILE STATISTICS:\n');
        fprintf('   EEG files per subject - Mean: %.1f, Range: %d-%d\n', ...
               mean(eegCounts), min(eegCounts), max(eegCounts));
        fprintf('   XDF files per subject - Mean: %.1f, Range: %d-%d\n', ...
               mean(xdfCounts), min(xdfCounts), max(xdfCounts));
        fprintf('   Pupil folders per subject - Mean: %.1f, Range: %d-%d\n', ...
               mean(pupilCounts), min(pupilCounts), max(pupilCounts));
    end
    
    fprintf('\n💾 Output directory: %s\n', outputPath);
    fprintf('==========================================================\n');
end

%% Helper Functions

function [matchedXdf, matchedCsv] = match_recording_events(xdfEvents, xdfTimes, csvEvents, csvTimes, verbose)
    % Match recording.begin and recording.end events
    
    matchedXdf = [];
    matchedCsv = [];
    
    if verbose
        fprintf('   🔧 Matching recording begin/end events...\n');
    end
    
    % Find begin and end events in both datasets
    recordingTypes = {'recording.begin', 'recording.end'};
    
    for i = 1:length(recordingTypes)
        eventType = recordingTypes{i};
        
        xdfIdx = find(strcmp(xdfEvents, eventType));
        csvIdx = find(strcmp(csvEvents, eventType));
        
        if ~isempty(xdfIdx) && ~isempty(csvIdx)
            % Use the first occurrence of each event type
            matchedXdf = [matchedXdf; xdfTimes(xdfIdx(1))];
            matchedCsv = [matchedCsv; csvTimes(csvIdx(1))];
            
            if verbose
                fprintf('      ✅ Matched %s: XDF=%.3fs, CSV=%.0fns\n', eventType, ...
                       xdfTimes(xdfIdx(1)), csvTimes(csvIdx(1)));
            end
        else
            if verbose
                fprintf('      ⚠️ %s not found in both datasets (XDF:%d, CSV:%d)\n', ...
                       eventType, length(xdfIdx), length(csvIdx));
            end
        end
    end
    
    if verbose
        fprintf('      ✅ Total matched recording events: %d\n', length(matchedXdf));
    end
end

function [matchedXdf, matchedCsv] = match_mixed_events(xdfEvents, xdfTimes, csvEvents, csvTimes, verbose)
    % Match mixed event types (both lsl.time_sync and recording events)
    
    if verbose
        fprintf('   🔧 Matching mixed event types...\n');
    end
    
    % First try LSL sync events
    [lslMatchedXdf, lslMatchedCsv] = match_sync_events_robust(xdfEvents, xdfTimes, csvEvents, csvTimes, false);
    
    % Then try recording events
    [recMatchedXdf, recMatchedCsv] = match_recording_events(xdfEvents, xdfTimes, csvEvents, csvTimes, false);
    
    % Combine results
    matchedXdf = [lslMatchedXdf; recMatchedXdf];
    matchedCsv = [lslMatchedCsv; recMatchedCsv];
    
    % Sort by time
    [matchedXdf, sortIdx] = sort(matchedXdf);
    matchedCsv = matchedCsv(sortIdx);
    
    if verbose
        fprintf('      ✅ Mixed matching: %d LSL sync + %d recording = %d total\n', ...
               length(lslMatchedXdf), length(recMatchedXdf), length(matchedXdf));
    end
end

function specialEvents = extract_special_events(eventsTable, verbose)
    % Extract special experimental events (not sync or recording events)
    
    % Define sync/recording events to exclude
    excludePatterns = {
        'recording.begin', 'recording.end', 'lsl.time_sync'
    };
    
    % Create mask for special events
    specialMask = true(height(eventsTable), 1);
    
    for i = 1:length(excludePatterns)
        pattern = excludePatterns{i};
        specialMask = specialMask & ~contains(eventsTable.name, pattern, 'IgnoreCase', true);
    end
    
    specialEvents = eventsTable(specialMask, :);
    
    if verbose && height(specialEvents) > 0
        fprintf('   🎯 Found %d special events:\n', height(specialEvents));
        
        % Show unique event types
        uniqueEvents = unique(specialEvents.name);
        for i = 1:min(10, length(uniqueEvents))  % Show first 10
            count = sum(strcmp(specialEvents.name, uniqueEvents{i}));
            fprintf('      - %s (%d)\n', uniqueEvents{i}, count);
        end
        if length(uniqueEvents) > 10
            fprintf('      ... and %d more types\n', length(uniqueEvents) - 10);
        end
    elseif verbose
        fprintf('   ⚠️ No special events found\n');
    end
end

function create_simple_events_timeline(eventsXdf, eegStart, eegEnd, verbose)
    % Create a simple timeline plot of special events
    
    if height(eventsXdf) == 0
        text(0.5, 0.5, 'No events to plot', 'Units', 'normalized', 'HorizontalAlignment', 'center');
        return;
    end
    
    % Get unique event types and assign colors
    uniqueTypes = unique(eventsXdf.name);
    colors = lines(length(uniqueTypes));
    
    hold on;
    
    % Plot each event type with a different color
    for i = 1:length(uniqueTypes)
        eventType = uniqueTypes{i};
        eventMask = strcmp(eventsXdf.name, eventType);
        eventTimes = eventsXdf.timestamp_xdf(eventMask);
        
        % Plot events as vertical lines
        for j = 1:length(eventTimes)
            line([eventTimes(j), eventTimes(j)], [i-0.4, i+0.4], ...
                 'Color', colors(i, :), 'LineWidth', 3);
        end
        
        % Add event type labels
        text(eegEnd + (eegEnd-eegStart)*0.01, i, strrep(eventType, '_', ' '), ...
             'FontSize', 8, 'Color', colors(i, :));
    end
    
    % Set axis properties
    xlim([eegStart, eegEnd + (eegEnd-eegStart)*0.15]);
    ylim([0.5, length(uniqueTypes) + 0.5]);
    set(gca, 'YTick', 1:length(uniqueTypes), 'YTickLabel', []);
    grid on;
    
    if verbose
        fprintf('      ✅ Plotted %d event types with %d total events\n', ...
               length(uniqueTypes), height(eventsXdf));
    end
end

function [eventCounts, eventTypes] = count_events_by_type(specialEvents)
    % Count events by type and sort by frequency
    
    [eventTypes, ~, idx] = unique(specialEvents.name);
    eventCounts = accumarray(idx, 1);
    
    % Sort by count (descending)
    [eventCounts, sortIdx] = sort(eventCounts, 'descend');
    eventTypes = eventTypes(sortIdx);
end

function create_detailed_special_events_figure(specialEvents, csvToXdfFun, eegStart, eegEnd, subjectName, outputPath, savePlots, verbose)
    % Create a detailed figure showing special events over time (simplified version)
    
    if height(specialEvents) == 0
        return;
    end
    
    try
        % Convert to XDF time
        specialEventsXdf = specialEvents;
        timestampCol = find_timestamp_column(specialEvents);
        specialEventsXdf.timestamp_xdf = csvToXdfFun(specialEvents.(timestampCol));
        
        % Filter to EEG time range
        timeFilter = specialEventsXdf.timestamp_xdf >= eegStart & specialEventsXdf.timestamp_xdf <= eegEnd;
        specialEventsFiltered = specialEventsXdf(timeFilter, :);
        
        if height(specialEventsFiltered) == 0
            if verbose
                fprintf('   ⚠️ No special events in EEG time range for detailed figure\n');
            end
            return;
        end
        
        fig2 = figure('Position', [150, 150, 1600, 800], 'Name', sprintf('%s Special Events Timeline', subjectName));
        
        % Main timeline plot
        subplot(2, 1, 1);
        create_simple_events_timeline(specialEventsFiltered, eegStart, eegEnd, false);
        title(sprintf('Special Events Timeline - %s (%d events)', subjectName, height(specialEventsFiltered)));
        xlabel('Time (seconds since LSL start)');
        
        % Event frequency over time
        subplot(2, 1, 2);
        if height(specialEventsFiltered) > 1
            % Create time bins (e.g., 1-minute bins)
            timeSpan = eegEnd - eegStart;
            binSize = max(60, timeSpan / 50);  % At least 60 seconds, max 50 bins
            timeBins = eegStart:binSize:eegEnd;
            
            eventCounts = histcounts(specialEventsFiltered.timestamp_xdf, timeBins);
            binCenters = timeBins(1:end-1) + binSize/2;
            
            bar(binCenters, eventCounts, 'FaceColor', [0.3, 0.6, 0.9], 'EdgeColor', 'k');
            xlabel('Time (seconds since LSL start)');
            ylabel('Events per bin');
            title(sprintf('Event Frequency (%.0fs bins)', binSize));
            grid on;
            xlim([eegStart, eegEnd]);
        else
            text(0.5, 0.5, 'Not enough events for frequency analysis', ...
                 'Units', 'normalized', 'HorizontalAlignment', 'center');
            title('Event Frequency');
        end
        
        if savePlots
            plotPath = fullfile(outputPath, sprintf('%s_special_events_detailed.png', subjectName));
            print(fig2, plotPath, '-dpng', '-r300');
            if verbose
                fprintf('   💾 Detailed special events plot saved: %s\n', plotPath);
            end
        end
        
        if verbose
            fprintf('   ✅ Detailed special events figure created (%d events)\n', height(specialEventsFiltered));
        end
        
    catch ME
        if verbose
            fprintf('   ❌ Error creating detailed special events figure: %s\n', ME.message);
        end
    end
end

function [matchedXdf, matchedCsv] = match_sync_events_robust(xdfEvents, xdfTimes, csvEvents, csvTimes, verbose)
    % Robust sync event matching using sequence numbers
    
    matchedXdf = [];
    matchedCsv = [];
    
    if verbose
        fprintf('   🔧 Matching sync events by sequence number...\n');
    end
    
    % Extract sequence numbers from event names
    xdfSeqNums = extract_sequence_numbers(xdfEvents);
    csvSeqNums = extract_sequence_numbers(csvEvents);
    
    if verbose
        fprintf('      XDF sequences: %s\n', mat2str(xdfSeqNums));
        fprintf('      CSV sequences: %s\n', mat2str(csvSeqNums));
    end
    
    % Find common sequence numbers
    [commonSeqs, xdfIdx, csvIdx] = intersect(xdfSeqNums, csvSeqNums);
    
    if length(commonSeqs) < 2
        if verbose
            fprintf('      ⚠️ Sequence matching failed, trying direct name matching...\n');
        end
        % Fallback to direct name matching
        [matchedXdf, matchedCsv] = match_sync_events_direct(xdfEvents, xdfTimes, csvEvents, csvTimes);
    else
        matchedXdf = xdfTimes(xdfIdx);
        matchedCsv = csvTimes(csvIdx);
        
        if verbose
            fprintf('      ✅ Matched %d events by sequence: %s\n', length(commonSeqs), mat2str(commonSeqs));
        end
    end
end

function seqNums = extract_sequence_numbers(eventNames)
    % Extract sequence numbers from LSL sync event names
    % e.g., 'lsl.time_sync.06a77f2e-3441-43d4-b310-4a3927636620.1' -> 1
    
    seqNums = [];
    
    for i = 1:length(eventNames)
        eventName = eventNames{i};
        
        % Look for pattern: lsl.time_sync.{uuid}.{number}
        tokens = regexp(eventName, 'lsl\.time_sync\.[^.]+\.(\d+)', 'tokens');
        
        if ~isempty(tokens) && ~isempty(tokens{1})
            seqNums(end+1) = str2double(tokens{1}{1});
        else
            % Try simpler pattern: just the last number after dots
            tokens = regexp(eventName, '\.(\d+)', 'tokens');
            if ~isempty(tokens) && ~isempty(tokens{1})
                seqNums(end+1) = str2double(tokens{1}{1});
            else
                % If no number found, use index as fallback
                seqNums(end+1) = i;
            end
        end
    end
end

function [matchedXdf, matchedCsv] = match_sync_events_direct(xdfEvents, xdfTimes, csvEvents, csvTimes)
    % Direct event name matching (fallback method)
    
    matchedXdf = [];
    matchedCsv = [];
    
    for i = 1:length(xdfEvents)
        xdfEvent = xdfEvents{i};
        for j = 1:length(csvEvents)
            csvEvent = csvEvents{j};
            if strcmp(xdfEvent, csvEvent)
                matchedXdf = [matchedXdf; xdfTimes(i)];
                matchedCsv = [matchedCsv; csvTimes(j)];
                break;
            end
        end
    end
end

function col_name = find_timestamp_column(table_data)
    % Find timestamp column in table with enhanced detection
    possible_names = {'timestampns', 'timestamp_ns_', 'timestamp [ns]', 'timestamp', 'startTimestampns', 'startTimestamp_ns_'};
    col_name = '';
    
    for i = 1:length(possible_names)
        % Try exact match first
        if ismember(possible_names{i}, table_data.Properties.VariableNames)
            col_name = possible_names{i};
            break;
        end
    end
    
    % If no exact match, try partial matching
    if isempty(col_name)
        for i = 1:length(possible_names)
            matches = contains(table_data.Properties.VariableNames, possible_names{i}, 'IgnoreCase', true);
            if any(matches)
                idx = find(matches, 1);
                col_name = table_data.Properties.VariableNames{idx};
                break;
            end
        end
    end
    
    % Final fallback: look for any column with 'timestamp' in the name
    if isempty(col_name)
        matches = contains(table_data.Properties.VariableNames, 'timestamp', 'IgnoreCase', true);
        if any(matches)
            idx = find(matches, 1);
            col_name = table_data.Properties.VariableNames{idx};
        end
    end
end

function col_name = find_column_by_names(col_names, possible_names)
    % Find column by possible names with enhanced flexibility
    col_name = '';
    
    % Try exact matches first
    for i = 1:length(possible_names)
        if any(strcmp(col_names, possible_names{i}))
            col_name = possible_names{i};
            break;
        end
    end
    
    % If no exact match, try case-insensitive partial matching
    if isempty(col_name)
        for i = 1:length(possible_names)
            matches = contains(col_names, possible_names{i}, 'IgnoreCase', true);
            if any(matches)
                idx = find(matches, 1);
                col_name = col_names{idx};
                break;
            end
        end
    end
    
    % Final fallback: try removing common suffixes/prefixes
    if isempty(col_name)
        for i = 1:length(possible_names)
            baseName = strrep(strrep(possible_names{i}, '_px', ''), 'px', '');
            matches = contains(col_names, baseName, 'IgnoreCase', true);
            if any(matches)
                idx = find(matches, 1);
                col_name = col_names{idx};
                break;
            end
        end
    end
end