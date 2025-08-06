%% EEG CONTINUOUS DATA ADDITION PIPELINE - Step 2
% Adds Faros ECG/ACC, Center of Mass, and Pupil Labs IMU data as continuous channels
% Run this after WaS4_00_pupilMerge.m
%
% This script:
% 1. Loads aligned EEG data from step 1
% 2. Loads Faros ECG/ACC data from EDF files
% 3. Loads Xsens Center of Mass data from MAT files
% 4. Loads Pupil Labs IMU data from IMU.csv files
% 5. Aligns all data to EEG sampling times using existing time mapping
% 6. Adds data as continuous channels to EEG structure
% 7. Saves enhanced EEG files

clear; clc; close all;

% Paths
alignedDataPath = '/Volumes/Work4TB/Seafile/WaS4/data/aligned/';  % Input: aligned data from step 1
outputPath = '/Volumes/Work4TB/Seafile/WaS4/data/continuous/';    % Output: EEG with continuous data
rawDataPath = '/Volumes/ergo/rohdaten/RohDaten_GRAIL/WaS4/DATA/';  % Raw data for Faros/Xsens

% Processing options
processAllSubjects = false;  % Set to false to process specific subjects
specificSubjects = [46];  % Only used if processAllSubjects = false
overwriteExisting = false;  % Set to true to reprocess existing files

% Interpolation and quality options
interpolationMethod = 'cubic';  % Method for interpolating continuous data
maxTimeGap = 1.0;              % Max gap (seconds) to interpolate across
minDataCoverage = 0.1;         % Minimum data coverage (10%) to add channel

% Output options
saveResults = true;
createPlots = true;
savePlots = true;
verboseOutput = true;

%% ==================== INITIALIZATION ====================

if verboseOutput
    fprintf('==========================================================\n');
    fprintf('EEG CONTINUOUS DATA ADDITION PIPELINE - Step 2\n');
    fprintf('==========================================================\n\n');
    fprintf('📁 Aligned data path: %s\n', alignedDataPath);
    fprintf('💾 Output path: %s\n', outputPath);
    fprintf('📂 Raw data path: %s\n\n', rawDataPath);
end

% Create output directory
if ~exist(outputPath, 'dir')
    mkdir(outputPath);
    if verboseOutput, fprintf('📁 Created output directory: %s\n\n', outputPath); end
end

% Find aligned EEG files
alignedFiles = dir(fullfile(alignedDataPath, '*_aligned.set'));
if isempty(alignedFiles)
    error('No aligned EEG files found in %s', alignedDataPath);
end

% Filter subjects if not processing all
if ~processAllSubjects
    validFiles = {};
    for s = 1:length(specificSubjects)
        subjectNum = specificSubjects(s);
        pattern = sprintf('WaS_%03d', subjectNum);
        for f = 1:length(alignedFiles)
            if contains(alignedFiles(f).name, pattern)
                validFiles{end+1} = alignedFiles(f).name;
                break;
            end
        end
    end
    if isempty(validFiles)
        error('No aligned files found for specified subjects');
    end
    subjectList = validFiles;
else
    subjectList = {alignedFiles.name};
end

if verboseOutput
    fprintf('🎯 Processing %d subjects:\n', length(subjectList));
    for i = 1:length(subjectList)
        fprintf('   - %s\n', strrep(subjectList{i}, '_aligned.set', ''));
    end
    fprintf('\n');
end

% Initialize tracking
processedSubjects = {};
failedSubjects = {};
skippedSubjects = {};
channelStats = [];

%% ==================== MAIN PROCESSING LOOP ====================

for s = 1:length(subjectList)
    alignedFile = subjectList{s};
    subjectFolder = strrep(alignedFile, '_aligned.set', '');
    
    if verboseOutput
        fprintf('==========================================================\n');
        fprintf('🔄 PROCESSING: %s (%d/%d)\n', subjectFolder, s, length(subjectList));
        fprintf('==========================================================\n\n');
    end
    
    try
        % Check if already processed
        outputFile = fullfile(outputPath, sprintf('%s_continuous.set', subjectFolder));
        if exist(outputFile, 'file') && ~overwriteExisting
            if verboseOutput
                fprintf('⏭️  Already processed (use overwriteExisting=true to reprocess)\n\n');
            end
            skippedSubjects{end+1} = subjectFolder;
            continue;
        end
        
        %% Step 1: Load aligned EEG data and time mapping
        if verboseOutput, fprintf('🧠 Loading aligned EEG data...\n'); end
        
        eegPath = fullfile(alignedDataPath, alignedFile);
        eeg = pop_loadset(eegPath);
        
        % Load alignment data for time mapping
        alignmentDataFile = fullfile(alignedDataPath, sprintf('%s_alignment_data.mat', subjectFolder));
        if ~exist(alignmentDataFile, 'file')
            error('Alignment data file not found: %s', alignmentDataFile);
        end
        
        alignmentData = load(alignmentDataFile);
        timeMapping = alignmentData.alignmentData.alignmentReport;
        
        % Create time conversion function (CSV nanoseconds to XDF seconds)
        csv_to_xdf_time = @(csv_time_ns) timeMapping.intercept + timeMapping.slope * (csv_time_ns * 1e-9);
        
        if verboseOutput
            fprintf('   ✅ EEG loaded: %d channels, %d samples, %.1f Hz\n', ...
                   eeg.nbchan, eeg.pnts, eeg.srate);
            fprintf('   ✅ Time mapping loaded: slope=%.9f, intercept=%.6f\n\n', ...
                   timeMapping.slope, timeMapping.intercept);
        end
        
        %% Step 2: Create EEG time vector in XDF time
        if verboseOutput, fprintf('⏱️  Creating EEG time vector...\n'); end
        
        % Find counter channel and create time vector
        cntIdx = find(strcmp({eeg.chanlocs.labels}, 'CNT'));
        if isempty(cntIdx)
            error('No CNT (counter) channel found in EEG data');
        end
        
        % Get original streams to find counter mapping
        streams = alignmentData.alignmentData.originalStreams;
        eegStreamIdx = [];
        for i = 1:length(streams)
            if contains(streams{i}.info.name, {'Counter', 'Sync', 'PROX'}, 'IgnoreCase', true)
                eegStreamIdx = i;
                break;
            end
        end
        
        if isempty(eegStreamIdx)
            error('Could not find EEG counter stream in alignment data');
        end
        
        % Match counter values between EEG and XDF to get time mapping
        [~, idxAmp, idxLSL] = intersect(eeg.data(cntIdx,:), streams{eegStreamIdx}.time_series);
        if length(idxAmp) < 10
            error('Too few matching counter values (%d) for reliable time mapping', length(idxAmp));
        end
        
        eegTimesXdf = NaN(1, eeg.pnts);
        eegTimesXdf(idxAmp) = streams{eegStreamIdx}.time_stamps(idxLSL);
        
        % Interpolate for missing time points
        validIdx = ~isnan(eegTimesXdf);
        eegTimesXdf = interp1(find(validIdx), eegTimesXdf(validIdx), 1:eeg.pnts, 'linear', 'extrap');
        
        eegTimeStart = min(eegTimesXdf);
        eegTimeEnd = max(eegTimesXdf);
        
        if verboseOutput
            fprintf('   ✅ EEG time mapping: %.3f - %.3f seconds (%.1f min)\n', ...
                   eegTimeStart, eegTimeEnd, (eegTimeEnd - eegTimeStart)/60);
            fprintf('   📊 Counter matches: %d/%d samples (%.1f%%)\n\n', ...
                   length(idxAmp), eeg.pnts, length(idxAmp)/eeg.pnts*100);
        end
        
        %% Step 3: Load and process Faros data
        if verboseOutput, fprintf('❤️  Loading Faros ECG/ACC data...\n'); end
        
        subjectRawPath = fullfile(rawDataPath, subjectFolder);
        farosFiles = dir(fullfile(subjectRawPath, '*.EDF'));
        
        channelsAdded = [];
        farosChannelsAdded = 0;
        
        if ~isempty(farosFiles)
            if verboseOutput, fprintf('   🔍 Found %d Faros files\n', length(farosFiles)); end
            
            farosDataMerged = [];
            
            % Process each Faros file
            for f = 1:length(farosFiles)
                filepath = fullfile(subjectRawPath, farosFiles(f).name);
                if verboseOutput, fprintf('   📥 Reading %s\n', farosFiles(f).name); end
                
                try
                    farosImport = edfread(filepath);
                    farosInfo = edfinfo(filepath);
                    
                    % Extract ECG and accelerometer data
                    if isfield(farosImport, 'ECG')
                        ecgData = cell2mat(farosImport.ECG);
                        ecgSampleRate = farosInfo.NumSamples(1) / farosInfo.DataRecordDuration;
                        ecgTime = (0:size(ecgData,1)-1) / ecgSampleRate;
                        
                        if verboseOutput, fprintf('      📊 ECG: %d samples at %.1f Hz\n', length(ecgData), ecgSampleRate); end
                    else
                        error('No ECG data found in Faros file');
                    end
                    
                    % Process accelerometer data
                    accData = [];
                    accLabels = {};
                    if isfield(farosImport, 'Accelerometer_X')
                        accX = cell2mat(farosImport.Accelerometer_X);
                        accY = cell2mat(farosImport.Accelerometer_Y);
                        accZ = cell2mat(farosImport.Accelerometer_Z);
                        accSampleRate = farosInfo.NumSamples(2) / farosInfo.DataRecordDuration;
                        
                        % Interpolate ACC to ECG sampling rate
                        accTime = (0:size(accX,1)-1) / accSampleRate;
                        accXInterp = interp1(accTime, accX, ecgTime, 'linear', 'extrap');
                        accYInterp = interp1(accTime, accY, ecgTime, 'linear', 'extrap');
                        accZInterp = interp1(accTime, accZ, ecgTime, 'linear', 'extrap');
                        
                        accData = [accXInterp', accYInterp', accZInterp'];
                        accLabels = {'FarosAccX', 'FarosAccY', 'FarosAccZ'};
                        
                        if verboseOutput, fprintf('      📊 ACC: %d samples at %.1f Hz (interpolated to ECG rate)\n', length(accX), accSampleRate); end
                    end
                    
                    % Create time-stamped data matrix
                    farosDataCurrent = [ecgTime', ecgData, accData];
                    
                    % Merge with existing data
                    if isempty(farosDataMerged)
                        farosDataMerged = farosDataCurrent;
                    else
                        % Simple concatenation - assumes files are in chronological order
                        timeOffset = farosDataMerged(end,1) + 1/ecgSampleRate;
                        farosDataCurrent(:,1) = farosDataCurrent(:,1) + timeOffset;
                        farosDataMerged = [farosDataMerged; farosDataCurrent];
                    end
                    
                catch ME
                    if verboseOutput
                        fprintf('      ❌ Error reading %s: %s\n', farosFiles(f).name, ME.message);
                    end
                    continue;
                end
            end
            
            % Process merged Faros data
            if ~isempty(farosDataMerged)
                % Convert relative time to XDF time using alignment
                % Use first available sync point from original LSL data
                farosLSLIdx = [];
                for i = 1:length(streams)
                    if contains(streams{i}.info.name, 'faros_ecg', 'IgnoreCase', true)
                        farosLSLIdx = i;
                        break;
                    end
                end
                
                if ~isempty(farosLSLIdx)
                    % Use cross-correlation to find time offset
                    farosLSLData = streams{farosLSLIdx}.time_series;
                    farosLSLTime = streams{farosLSLIdx}.time_stamps;
                    
                    % Resample LSL data to match Faros sampling rate
                    lslECGResampled = interp1(farosLSLTime, farosLSLData, ...
                                            farosLSLTime(1):1/ecgSampleRate:farosLSLTime(end), 'pchip');
                    
                    % Find best alignment using cross-correlation
                    [r, lags] = xcorr(farosDataMerged(:,2), lslECGResampled);
                    [~, maxIdx] = max(r);
                    offsetSamples = lags(maxIdx);
                    offsetTime = offsetSamples / ecgSampleRate;
                    
                    % Apply time offset
                    farosDataMerged(:,1) = farosDataMerged(:,1) + farosLSLTime(1) + offsetTime;
                    
                    if verboseOutput
                        fprintf('   🔧 Faros time alignment: offset = %.3f seconds (r_max = %.3f)\n', offsetTime, max(r));
                    end
                else
                    % Fallback: use EEG start time
                    farosDataMerged(:,1) = farosDataMerged(:,1) + eegTimeStart;
                    if verboseOutput, fprintf('   ⚠️ Using fallback time alignment (no LSL Faros reference)\n'); end
                end
                
                % Filter to EEG time range
                timeFilter = farosDataMerged(:,1) >= eegTimeStart & farosDataMerged(:,1) <= eegTimeEnd;
                farosFiltered = farosDataMerged(timeFilter, :);
                
                if ~isempty(farosFiltered)
                    % Interpolate to EEG sample times
                    ecgInterp = interp1(farosFiltered(:,1), farosFiltered(:,2), eegTimesXdf, interpolationMethod, NaN);
                    
                    % Check coverage
                    ecgCoverage = sum(~isnan(ecgInterp)) / length(ecgInterp);
                    
                    if ecgCoverage >= minDataCoverage
                        % Add ECG channel
                        newIdx = eeg.nbchan + 1;
                        eeg.data(newIdx, :) = ecgInterp;
                        eeg.chanlocs(newIdx).labels = 'FarosECG';
                        eeg.chanlocs(newIdx).type = 'ECG';
                        eeg.nbchan = newIdx;
                        channelsAdded(end+1) = newIdx;
                        farosChannelsAdded = farosChannelsAdded + 1;
                        
                        if verboseOutput
                            fprintf('   ✅ Added ECG channel (%.1f%% coverage)\n', ecgCoverage * 100);
                        end
                    else
                        if verboseOutput
                            fprintf('   ⚠️ ECG coverage too low (%.1f%% < %.1f%%)\n', ecgCoverage * 100, minDataCoverage * 100);
                        end
                    end
                    
                    % Add accelerometer channels
                    if size(farosFiltered, 2) > 2
                        for accCh = 1:3
                            if size(farosFiltered, 2) > 2 + accCh - 1
                                accInterp = interp1(farosFiltered(:,1), farosFiltered(:,2+accCh), eegTimesXdf, interpolationMethod, NaN);
                                accCoverage = sum(~isnan(accInterp)) / length(accInterp);
                                
                                if accCoverage >= minDataCoverage
                                    newIdx = eeg.nbchan + 1;
                                    eeg.data(newIdx, :) = accInterp;
                                    eeg.chanlocs(newIdx).labels = accLabels{accCh};
                                    eeg.chanlocs(newIdx).type = 'ACC';
                                    eeg.nbchan = newIdx;
                                    channelsAdded(end+1) = newIdx;
                                    farosChannelsAdded = farosChannelsAdded + 1;
                                    
                                    if verboseOutput
                                        fprintf('   ✅ Added %s channel (%.1f%% coverage)\n', accLabels{accCh}, accCoverage * 100);
                                    end
                                else
                                    if verboseOutput
                                        fprintf('   ⚠️ %s coverage too low (%.1f%% < %.1f%%)\n', accLabels{accCh}, accCoverage * 100, minDataCoverage * 100);
                                    end
                                end
                            end
                        end
                    end
                else
                    if verboseOutput, fprintf('   ⚠️ No Faros data in EEG time range\n'); end
                end
            else
                if verboseOutput, fprintf('   ⚠️ No valid Faros data loaded\n'); end
            end
        else
            if verboseOutput, fprintf('   ⚠️ No Faros files found\n'); end
        end
        
        if verboseOutput, fprintf('   📊 Added %d Faros channels\n\n', farosChannelsAdded); end
        
        %% Step 4: Load and process Xsens Center of Mass data
        if verboseOutput, fprintf('🤸 Loading Xsens Center of Mass data...\n'); end
        
        xsensFiles = dir(fullfile(subjectRawPath, '*xsens*.mat'));
        xsensChannelsAdded = 0;
        
        if ~isempty(xsensFiles)
            if verboseOutput, fprintf('   🔍 Found %d Xsens files\n', length(xsensFiles)); end
            
            try
                xsensPath = fullfile(subjectRawPath, xsensFiles(1).name);
                if verboseOutput, fprintf('   📥 Reading %s\n', xsensFiles(1).name); end
                
                xsensData = load(xsensPath);
                
                if isfield(xsensData, 'xsensData') && isfield(xsensData.xsensData, 'center_of_mass')
                    xsensCoMData = [];
                    xsensTimeData = [];
                    
                    % Process each recording segment
                    for x = 1:length(xsensData.xsensData.center_of_mass)
                        frameRate = xsensData.xsensData.info{x}{5,2};
                        comData = xsensData.xsensData.center_of_mass{x};
                        
                        % Create time vector
                        segmentTime = (0:height(comData)-1) / frameRate;
                        
                        % Convert to XDF time using LSL reference
                        xsensLSLIdx = [];
                        for i = 1:length(streams)
                            if contains(streams{i}.info.name, 'CenterOfMass', 'IgnoreCase', true)
                                xsensLSLIdx = i;
                                break;
                            end
                        end
                        
                        if ~isempty(xsensLSLIdx)
                            % Use cross-correlation for time alignment
                            lslCoMX = streams{xsensLSLIdx}.time_series(1,:);
                            lslTime = streams{xsensLSLIdx}.time_stamps;
                            
                            % Resample LSL data
                            lslCoMResampled = interp1(lslTime, lslCoMX, lslTime(1):1/frameRate:lslTime(end), 'pchip');
                            
                            % Find alignment
                            [r, lags] = xcorr(lslCoMResampled, comData.CoMPosX);
                            [~, maxIdx] = max(r);
                            offsetTime = lags(maxIdx) / frameRate;
                            
                            segmentTime = segmentTime + lslTime(1) + offsetTime;
                            
                            if verboseOutput
                                fprintf('      🔧 Segment %d alignment: offset = %.3f s, r_max = %.3f\n', x, offsetTime, max(r));
                            end
                        else
                            % Fallback alignment
                            segmentTime = segmentTime + eegTimeStart + x * 60; % Assume 1-minute segments
                            if verboseOutput, fprintf('      ⚠️ Using fallback alignment for segment %d\n', x); end
                        end
                        
                        % Add to merged data
                        xsensTimeData = [xsensTimeData, segmentTime];
                        
                        comArray = table2array(comData(:, 2:end)); % Skip time column
                        if isempty(xsensCoMData)
                            xsensCoMData = comArray';
                        else
                            xsensCoMData = [xsensCoMData, comArray'];
                        end
                    end
                    
                    % Filter to EEG time range
                    timeFilter = xsensTimeData >= eegTimeStart & xsensTimeData <= eegTimeEnd;
                    xsensTimeFiltered = xsensTimeData(timeFilter);
                    xsensDataFiltered = xsensCoMData(:, timeFilter);
                    
                    if ~isempty(xsensTimeFiltered)
                        % Get channel names from the first segment
                        comChannelNames = xsensData.xsensData.center_of_mass{1}.Properties.VariableNames(2:end);
                        
                        % Add each CoM channel
                        for ch = 1:size(xsensDataFiltered, 1)
                            comInterp = interp1(xsensTimeFiltered, xsensDataFiltered(ch, :), eegTimesXdf, interpolationMethod, NaN);
                            comCoverage = sum(~isnan(comInterp)) / length(comInterp);
                            
                            if comCoverage >= minDataCoverage
                                newIdx = eeg.nbchan + 1;
                                eeg.data(newIdx, :) = comInterp;
                                eeg.chanlocs(newIdx).labels = ['CoM_' comChannelNames{ch}];
                                eeg.chanlocs(newIdx).type = 'Xsens';
                                eeg.nbchan = newIdx;
                                channelsAdded(end+1) = newIdx;
                                xsensChannelsAdded = xsensChannelsAdded + 1;
                                
                                if verboseOutput
                                    fprintf('   ✅ Added %s channel (%.1f%% coverage)\n', comChannelNames{ch}, comCoverage * 100);
                                end
                            else
                                if verboseOutput
                                    fprintf('   ⚠️ %s coverage too low (%.1f%% < %.1f%%)\n', comChannelNames{ch}, comCoverage * 100, minDataCoverage * 100);
                                end
                            end
                        end
                    else
                        if verboseOutput, fprintf('   ⚠️ No Xsens data in EEG time range\n'); end
                    end
                else
                    if verboseOutput, fprintf('   ⚠️ Invalid Xsens data structure\n'); end
                end
                
            catch ME
                if verboseOutput
                    fprintf('   ❌ Error loading Xsens data: %s\n', ME.message);
                end
            end
        else
            if verboseOutput, fprintf('   ⚠️ No Xsens files found\n'); end
        end
        
        if verboseOutput, fprintf('   📊 Added %d Xsens CoM channels\n\n', xsensChannelsAdded); end
        
        %% Step 5: Load and process Pupil Labs IMU data
        if verboseOutput, fprintf('📱 Loading Pupil Labs IMU data...\n'); end
        
        pupilDirs = dir(fullfile(subjectRawPath, '*pupil*'));
        pupilDirs = pupilDirs([pupilDirs.isdir]);
        imuChannelsAdded = 0;
        
        if ~isempty(pupilDirs)
            if verboseOutput, fprintf('   🔍 Found %d Pupil Labs folders\n', length(pupilDirs)); end
            
            pupilIMUData = [];
            
            % Process all pupil directories
            for pupilIdx = 1:length(pupilDirs)
                pupilDir = fullfile(subjectRawPath, pupilDirs(pupilIdx).name);
                if verboseOutput, fprintf('   📂 Processing folder: %s\n', pupilDirs(pupilIdx).name); end
                
                % Find actual data folder
                dataFolders = dir(pupilDir);
                dataFolders = dataFolders([dataFolders.isdir] & ~strcmp({dataFolders.name}, '.') & ~strcmp({dataFolders.name}, '..'));
                
                if ~isempty(dataFolders)
                    [~, latest_idx] = max([dataFolders.datenum]);
                    pupilDataPath = fullfile(pupilDir, dataFolders(latest_idx).name);
                else
                    pupilDataPath = pupilDir;
                end
                
                % Look for IMU.csv
                imuPath = fullfile(pupilDataPath, 'IMU.csv');
                if exist(imuPath, 'file')
                    try
                        if verboseOutput, fprintf('      📥 Reading IMU.csv\n'); end
                        
                        imuData = readtable(imuPath);
                        
                        if height(imuData) > 0
                            if verboseOutput, fprintf('         📊 IMU data: %d samples\n', height(imuData)); end
                            
                            % Find timestamp column
                            timestampCol = find_timestamp_column_imu(imuData);
                            if isempty(timestampCol)
                                if verboseOutput, fprintf('         ⚠️ No timestamp column found\n'); end
                                continue;
                            end
                            
                            % Convert timestamps to XDF time
                            imuData.timestamp_xdf = csv_to_xdf_time(imuData.(timestampCol));
                            
                            % Merge with existing data
                            if isempty(pupilIMUData)
                                pupilIMUData = imuData;
                            else
                                pupilIMUData = vertcat(pupilIMUData, imuData);
                            end
                        else
                            if verboseOutput, fprintf('         ⚠️ Empty IMU data\n'); end
                        end
                        
                    catch ME
                        if verboseOutput
                            fprintf('         ❌ Error reading IMU.csv: %s\n', ME.message);
                        end
                    end
                else
                    if verboseOutput, fprintf('      ⚠️ IMU.csv not found\n'); end
                end
            end
            
            % Process merged IMU data
            if ~isempty(pupilIMUData)
                % Sort by time
                pupilIMUData = sortrows(pupilIMUData, 'timestamp_xdf');
                
                % Filter to EEG time range
                timeFilter = pupilIMUData.timestamp_xdf >= eegTimeStart & pupilIMUData.timestamp_xdf <= eegTimeEnd;
                imuFiltered = pupilIMUData(timeFilter, :);
                
                if height(imuFiltered) > 0
                    % Define IMU channels to extract
                    imuChannels = {
                        'gyro_x_rad_s_', 'IMU_GyroX';
                        'gyro_y_rad_s_', 'IMU_GyroY';
                        'gyro_z_rad_s_', 'IMU_GyroZ';
                        'accel_x_m_s2_', 'IMU_AccX';
                        'accel_y_m_s2_', 'IMU_AccY';
                        'accel_z_m_s2_', 'IMU_AccZ';
                        'mag_x_uT_', 'IMU_MagX';
                        'mag_y_uT_', 'IMU_MagY';
                        'mag_z_uT_', 'IMU_MagZ'
                    };
                    
                    % Add each IMU channel
                    for ch = 1:size(imuChannels, 1)
                        csvColName = imuChannels{ch, 1};
                        eegChLabel = imuChannels{ch, 2};
                        
                        % Try flexible column name matching
                        actualCol = find_imu_column(imuFiltered.Properties.VariableNames, csvColName);
                        
                        if ~isempty(actualCol)
                            imuInterp = interp1(imuFiltered.timestamp_xdf, imuFiltered.(actualCol), eegTimesXdf, interpolationMethod, NaN);
                            imuCoverage = sum(~isnan(imuInterp)) / length(imuInterp);
                            
                            if imuCoverage >= minDataCoverage
                                newIdx = eeg.nbchan + 1;
                                eeg.data(newIdx, :) = imuInterp;
                                eeg.chanlocs(newIdx).labels = eegChLabel;
                                eeg.chanlocs(newIdx).type = 'IMU';
                                eeg.nbchan = newIdx;
                                channelsAdded(end+1) = newIdx;
                                imuChannelsAdded = imuChannelsAdded + 1;
                                
                                if verboseOutput
                                    fprintf('   ✅ Added %s channel (%.1f%% coverage)\n', eegChLabel, imuCoverage * 100);
                                end
                            else
                                if verboseOutput
                                    fprintf('   ⚠️ %s coverage too low (%.1f%% < %.1f%%)\n', eegChLabel, imuCoverage * 100, minDataCoverage * 100);
                                end
                            end
                        else
                            if verboseOutput
                                fprintf('   ⚠️ IMU column %s not found\n', csvColName);
                            end
                        end
                    end
                else
                    if verboseOutput, fprintf('   ⚠️ No IMU data in EEG time range\n'); end
                end
            else
                if verboseOutput, fprintf('   ⚠️ No IMU data loaded\n'); end
            end
        else
            if verboseOutput, fprintf('   ⚠️ No Pupil Labs folders found\n'); end
        end
        
        if verboseOutput, fprintf('   📊 Added %d IMU channels\n\n', imuChannelsAdded); end
        
        %% Step 6: Create quality plots
        if createPlots
            if verboseOutput, fprintf('📈 Creating quality plots...\n'); end
            
            fig = figure('Position', [100, 100, 1600, 1000], 'Name', sprintf('%s - Continuous Data Quality', subjectFolder));
            
            totalNewChannels = farosChannelsAdded + xsensChannelsAdded + imuChannelsAdded;
            
            if totalNewChannels > 0
                % Plot 1: Data coverage for new channels
                subplot(2, 3, 1);
                if ~isempty(channelsAdded)
                    coverage = [];
                    labels = {};
                    for ch = channelsAdded
                        validData = ~isnan(eeg.data(ch, :));
                        coverage(end+1) = sum(validData) / length(validData) * 100;
                        labels{end+1} = eeg.chanlocs(ch).labels;
                    end
                    
                    barh(coverage);
                    set(gca, 'YTickLabel', labels, 'YTick', 1:length(labels));
                    xlabel('Coverage (%)');
                    title('Data Coverage by Channel');
                    grid on;
                end
                
                % Plot 2: Sample data from new channels
                subplot(2, 3, 2);
                if ~isempty(channelsAdded) && length(channelsAdded) >= 1
                    % Show first new channel
                    ch = channelsAdded(1);
                    timeAxis = (1:eeg.pnts) / eeg.srate;
                    plot(timeAxis, eeg.data(ch, :));
                    xlabel('Time (s)');
                    ylabel('Amplitude');
                    title(sprintf('Sample Data: %s', eeg.chanlocs(ch).labels));
                    grid on;
                end
                
                % Plot 3: Channel type summary
                subplot(2, 3, 3);
                channelTypes = {};
                for ch = channelsAdded
                    channelTypes{end+1} = eeg.chanlocs(ch).type;
                end
                if ~isempty(channelTypes)
                    [uniqueTypes, ~, idx] = unique(channelTypes);
                    counts = accumarray(idx, 1);
                    pie(counts, uniqueTypes);
                    title('Added Channels by Type');
                end
                
                % Plot 4-6: Time series examples for each data type
                plotIdx = 4;
                dataTypes = {'ECG', 'Xsens', 'IMU'};
                for dt = 1:length(dataTypes)
                    subplot(2, 3, plotIdx);
                    plotIdx = plotIdx + 1;
                    
                    % Find channels of this type
                    typeChannels = [];
                    for ch = channelsAdded
                        if strcmp(eeg.chanlocs(ch).type, dataTypes{dt}) || contains(eeg.chanlocs(ch).labels, dataTypes{dt})
                            typeChannels(end+1) = ch;
                        end
                    end
                    
                    if ~isempty(typeChannels)
                        timeAxis = (1:eeg.pnts) / eeg.srate;
                        hold on;
                        colors = lines(length(typeChannels));
                        for i = 1:min(3, length(typeChannels)) % Show max 3 channels
                            ch = typeChannels(i);
                            plot(timeAxis, eeg.data(ch, :), 'Color', colors(i, :), 'DisplayName', eeg.chanlocs(ch).labels);
                        end
                        xlabel('Time (s)');
                        ylabel('Amplitude');
                        title(sprintf('%s Data', dataTypes{dt}));
                        legend('Location', 'best');
                        grid on;
                    else
                        text(0.5, 0.5, sprintf('No %s data added', dataTypes{dt}), ...
                             'Units', 'normalized', 'HorizontalAlignment', 'center');
                        title(sprintf('%s Data', dataTypes{dt}));
                    end
                end
            else
                % No channels added
                text(0.5, 0.5, 'No continuous channels were added', ...
                     'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 16);
                title('No Data Added');
            end
            
            if savePlots && totalNewChannels > 0
                plotPath = fullfile(outputPath, sprintf('%s_continuous_quality.png', subjectFolder));
                print(fig, plotPath, '-dpng', '-r300');
                if verboseOutput, fprintf('   💾 Quality plot saved: %s\n', plotPath); end
            end
            
            if verboseOutput, fprintf('   ✅ Quality plots created\n\n'); end
        end
        
        %% Step 7: Update EEG structure and save
        if verboseOutput, fprintf('💾 Updating EEG structure and saving...\n'); end
        
        % Update EEG structure
        eeg = eeg_checkset(eeg);
        
        % Save enhanced EEG file
        if saveResults
            pop_saveset(eeg, 'filename', sprintf('%s_continuous.set', subjectFolder), 'filepath', outputPath);
            
            % Save processing info
            continuousInfo = struct();
            continuousInfo.originalChannels = eeg.nbchan - length(channelsAdded);
            continuousInfo.addedChannels = length(channelsAdded);
            continuousInfo.farosChannels = farosChannelsAdded;
            continuousInfo.xsensChannels = xsensChannelsAdded;
            continuousInfo.imuChannels = imuChannelsAdded;
            continuousInfo.addedChannelIndices = channelsAdded;
            continuousInfo.processingOptions.interpolationMethod = interpolationMethod;
            continuousInfo.processingOptions.maxTimeGap = maxTimeGap;
            continuousInfo.processingOptions.minDataCoverage = minDataCoverage;
            continuousInfo.processingInfo.timestamp = datetime('now');
            continuousInfo.processingInfo.matlabVersion = version;
            continuousInfo.processingInfo.script = mfilename;
            
            infoPath = fullfile(outputPath, sprintf('%s_continuous_info.mat', subjectFolder));
            save(infoPath, 'continuousInfo', '-v7.3');
            
            if verboseOutput
                fprintf('   ✅ EEG data saved: %s\n', outputFile);
                fprintf('   ✅ Processing info saved: %s\n', infoPath);
            end
        end
        
        % Record statistics
        channelStats = [channelStats; struct('subject', subjectFolder, ...
                                            'originalChannels', eeg.nbchan - length(channelsAdded), ...
                                            'addedChannels', length(channelsAdded), ...
                                            'farosChannels', farosChannelsAdded, ...
                                            'xsensChannels', xsensChannelsAdded, ...
                                            'imuChannels', imuChannelsAdded)];
        
        processedSubjects{end+1} = subjectFolder;
        
        if verboseOutput
            fprintf('\n✅ %s COMPLETED SUCCESSFULLY!\n', subjectFolder);
            fprintf('   Added channels: %d total (Faros:%d, Xsens:%d, IMU:%d)\n', ...
                   length(channelsAdded), farosChannelsAdded, xsensChannelsAdded, imuChannelsAdded);
            fprintf('   Final EEG structure: %d channels, %d samples\n', eeg.nbchan, eeg.pnts);
        end
        
    catch ME
        failedSubjects{end+1} = sprintf('%s: %s', subjectFolder, ME.message);
        
        if verboseOutput
            fprintf('\n❌ ERROR processing %s:\n', subjectFolder);
            fprintf('   %s\n', ME.message);
            if length(ME.stack) > 0
                fprintf('   Line %d in %s\n', ME.stack(1).line, ME.stack(1).name);
            end
        end
        continue;
    end
    
    % Close figures
    close all;
end

%% ==================== FINAL SUMMARY ====================

if verboseOutput
    fprintf('\n==========================================================\n');
    fprintf('CONTINUOUS DATA ADDITION PIPELINE SUMMARY\n');
    fprintf('==========================================================\n\n');
    
    fprintf('✅ SUCCESSFULLY PROCESSED (%d subjects):\n', length(processedSubjects));
    for i = 1:length(processedSubjects)
        if i <= length(channelStats)
            stats = channelStats(i);
            fprintf('   %s: %d→%d channels (+%d: Faros:%d, Xsens:%d, IMU:%d)\n', ...
                   stats.subject, stats.originalChannels, ...
                   stats.originalChannels + stats.addedChannels, stats.addedChannels, ...
                   stats.farosChannels, stats.xsensChannels, stats.imuChannels);
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
    if ~isempty(channelStats)
        allAdded = [channelStats.addedChannels];
        allFaros = [channelStats.farosChannels];
        allXsens = [channelStats.xsensChannels];
        allIMU = [channelStats.imuChannels];
        
        fprintf('\n📊 OVERALL STATISTICS:\n');
        fprintf('   Channels added per subject - Mean: %.1f, Range: %d-%d\n', ...
               mean(allAdded), min(allAdded), max(allAdded));
        fprintf('   Faros channels - Mean: %.1f, Total: %d\n', mean(allFaros), sum(allFaros));
        fprintf('   Xsens channels - Mean: %.1f, Total: %d\n', mean(allXsens), sum(allXsens));
        fprintf('   IMU channels - Mean: %.1f, Total: %d\n', mean(allIMU), sum(allIMU));
    end
    
    fprintf('\n💾 Output directory: %s\n', outputPath);
    fprintf('==========================================================\n');
end

%% ==================== HELPER FUNCTIONS ====================

function timestampCol = find_timestamp_column_imu(tableData)
    % Find timestamp column in IMU data
    possibleNames = {'timestamp_ns_', 'timestampns', 'timestamp [ns]', 'timestamp'};
    timestampCol = '';
    
    for i = 1:length(possibleNames)
        if ismember(possibleNames{i}, tableData.Properties.VariableNames)
            timestampCol = possibleNames{i};
            break;
        end
    end
    
    % Try partial matching
    if isempty(timestampCol)
        for i = 1:length(possibleNames)
            matches = contains(tableData.Properties.VariableNames, possibleNames{i}, 'IgnoreCase', true);
            if any(matches)
                idx = find(matches, 1);
                timestampCol = tableData.Properties.VariableNames{idx};
                break;
            end
        end
    end
end

function colName = find_imu_column(colNames, targetName)
    % Find IMU column with flexible matching
    colName = '';
    
    % Try exact match first
    if any(strcmp(colNames, targetName))
        colName = targetName;
        return;
    end
    
    % Try partial matching
    matches = contains(colNames, targetName, 'IgnoreCase', true);
    if any(matches)
        idx = find(matches, 1);
        colName = colNames{idx};
        return;
    end
    
    % Try removing underscores and matching
    targetClean = strrep(targetName, '_', '');
    for i = 1:length(colNames)
        colClean = strrep(colNames{i}, '_', '');
        if contains(colClean, targetClean, 'IgnoreCase', true)
            colName = colNames{i};
            break;
        end
    end
end