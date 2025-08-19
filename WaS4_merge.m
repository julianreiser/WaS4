function [eegMerged, streamsMerged] = WaS4_merge(subject, dataPath, outPath)

folders = dir(fullfile(dataPath, ['*WaS*' sprintf('%03d', subject) '*']));

if isempty(folders)
    warning(['no data found for subject ' sprintf('%03d', subject)]);
    return
end

subjectFolder = folders(1).name;

%% gather all required data

%% load eeg data from bdf file
files = dir(fullfile(dataPath, subjectFolder,  '*.bdf'));
ampTmp = {};

if isempty(files)
    warning('No EEG files found, stopping...')
    return
end
for f = 1:length(files)
    filepath = fullfile(dataPath, subjectFolder, files(f).name);
    disp(['reading ' filepath '...'])
    ampTmp{f} = pop_biosig(filepath);
end

% merge files if > 1
amp = ampTmp{1};
cntIdx = find(strcmp({amp.chanlocs.labels}, 'CNT'));
ampCNTOffsets = [0];
for f = 2:length(ampTmp)
    newTime = amp.times(end) + ampTmp{f}.times + (1000/amp.srate);
    amp.times = [amp.times newTime];

    ampTmp{f}.data(cntIdx, :) = ampTmp{f}.data(cntIdx, :) + amp.data(cntIdx,end); % adjust counter value to continue from the last file
    ampCNTOffsets = [ampCNTOffsets amp.data(cntIdx, end)]; % offsets to adjust counter values in LSL data later (last counter value before new recording)

    amp.data = [amp.data, ampTmp{f}.data];
    amp.pnts = amp.pnts + ampTmp{f}.pnts;
end
disp('--- done reading EEG ---')
% amp <- eeg data


%% load faros data from edf file
files = dir(fullfile(dataPath, subjectFolder, '*.EDF'));
farosData = {};

if isempty(files)
    warning('No Faros files found, stopping...')
    return
end

% read files
for f = 1:length(files)
    filepath = fullfile(dataPath, subjectFolder, files(f).name);
    disp(['reading ' filepath '...'])
    farosImport = edfread(filepath);
    farosInfo = edfinfo(filepath);
    ecgSampRate = 1/farosInfo.NumSamples(1);

    newECGData = cell2mat(farosImport.ECG);
    newACCData = [cell2mat(farosImport.Accelerometer_X), cell2mat(farosImport.Accelerometer_Y), cell2mat(farosImport.Accelerometer_Z)];
    newECGTime = (0:size(newECGData, 1)-1) * (1/farosInfo.NumSamples(1));
    newACCTime = (0:size(newACCData, 1)-1) * (1/farosInfo.NumSamples(2));

    newFarosData = [newECGTime', newECGData, zeros(length(newECGTime),size(newACCData,2))];
    for a = 1:size(newACCData,2)
        newFarosData(:, 2+a) = interp1(newACCTime, newACCData(:, a), newECGTime, 'makima');
    end

    farosData{f} = newFarosData;
end
disp('--- done reading ECG ---')
% farosData <- faros timestamps, eeg, acc X/Y/Z


%% load lsl data from xdf file
files = dir(fullfile(dataPath, subjectFolder, '*.xdf'));
streamsTmp = {};

if isempty(files)
    warning('No LSL data found, stopping...');
    return
end

for f = 1:length(files)
    filepath = fullfile(dataPath, subjectFolder, files(f).name);
    disp(['reading ' filepath '...'])
    streamsTmp{f} = load_xdf(filepath);
end

% if files > 1 look for timestamp reset (labrecorder pc restart)
if length(streamsTmp) > 1
    for f = 2:length(streamsTmp)
        reset = false;

        idx1 = cellfun(@(x) contains(x.info.name, {'Counter', 'Sync', 'PROX'}, 'IgnoreCase', true), streamsTmp{f-1});
        idx2 = cellfun(@(x) contains(x.info.name, {'Counter', 'Sync', 'PROX'}, 'IgnoreCase', true), streamsTmp{f});

        if streamsTmp{f}{idx2}.time_stamps(1) < streamsTmp{f-1}{idx1}.time_stamps(end)
            reset = true;
        end

        % fix timestamps
        if reset
            warning('Detected timestamp reset in LSL data, attempting to fix...');
            % use eeg counter if possible, ecg otherwise
            if length(ampTmp) == 1
                idxLast = cellfun(@(x) contains(x.info.name, {'Counter', 'Sync', 'PROX'}, 'IgnoreCase', true), streamsTmp{f-1});
                idxCurrent = cellfun(@(x) contains(x.info.name, {'Counter', 'Sync', 'PROX'}, 'IgnoreCase', true), streamsTmp{f});

                lastEnd = amp.time(amp.data(cntIdx,:) == streamsTmp{f-1}{idxLast}.time_series(end));
                currentFirst = amp.time(amp.data(cntIdx,:) == streamsTmp{f}{idxCurrent}.time_series(1));

                lslTimeOffset = streamsTmp{f-1}{idxLast}.time_stamps(end) + (currentFirst - lastEnd)/1000 - streamsTmp{f}{idxCurrent}.time_stamps(1);
            elseif length(farosData) == 1
                idxLast = cellfun(@(x) contains(x.info.name, 'faros_ecg'), streamsTmp{f-1});
                idxCurrent = cellfun(@(x) contains(x.info.name, 'faros_ecg'), streamsTmp{f});

                newLSLTimeLast = streamsTmp{f-1}{idxLast}.time_stamps(1):ecgSampRate:streamsTmp{f-1}{idxLast}.time_stamps(end);
                newLSLECGLast = interp1(streamsTmp{f-1}{idxLast}.time_stamps, double(streamsTmp{f-1}{idxLast}.time_series), newLSLTimeLast, 'pchip');
                [r,lags] = xcorr(farosData{1}(:,2), newLSLECGLast);
                [~, maxIdx] = max(r);
                offsetLast = lags(maxIdx);

                newLSLTimeCurrent = streamsTmp{f}{idxCurrent}.time_stamps(1):ecgSampRate:streamsTmp{f}{idxCurrent}.time_stamps(end);
                newLSLECGCurrent = interp1(streamsTmp{f}{idxCurrent}.time_stamps, double(streamsTmp{f}{idxCurrent}.time_series), newLSLTimeCurrent, 'pchip');
                [r,lags] = xcorr(farosData{1}(:,2), newLSLECGCurrent);
                [~, maxIdx] = max(r);
                offsetCurrent = lags(maxIdx);

                lslTimeOffset = streamsTmp{f-1}{idxLast}.time_stamps(1) + (offsetCurrent - offsetLast)*ecgSampRate - streamsTmp{f}{idxCurrent}.time_stamps(1);
            else
                warning("Could not fix LSL timestamps, stopping...");
                return
            end

            % use offset to adjust times for all parts of LSL data
            for nf = f:length(streamsTmp)
                for d = 1:length(streamsTmp{nf})
                    if ~isempty(streamsTmp{nf}{d}.time_stamps)
                        streamsTmp{nf}{d}.time_stamps = streamsTmp{nf}{d}.time_stamps + lslTimeOffset;
                    end
                end
            end
        end
    end

end

% merge files if > 1
streams = streamsTmp{1};
for f = 2:length(streamsTmp)
    nextStream = streamsTmp{f};
    for n = 1:length(nextStream)
        for id = 1:length(streams)
            if strcmp(streams{id}.info.name, nextStream{n}.info.name)
                if strcmp(streams{id}.info.type, "EEG") && nextStream{n}.time_series(1) < streams{id}.time_series(end)
                    % adjust EEG counter values
                    offsetIdx = find(ampCNTOffsets + nextStream{n}.time_series(1) > streams{id}.time_series(end), 1, 'first');
                    nextStream{n}.time_series = nextStream{n}.time_series + ampCNTOffsets(offsetIdx);
                end
                streams{id}.time_stamps = [streams{id}.time_stamps nextStream{n}.time_stamps];
                streams{id}.time_series = [streams{id}.time_series nextStream{n}.time_series];
            end
        end
    end
end
% streamNames = {};
% for id = 1:length(streams)
%     streamNames{id} = streams{id}.info.name;
% end
disp('--- done reading LSL ---')
% streams <- lsl recorder data

%% load xsens data from mat file
% run convertXsens first to generate .mat
xsensData = [];
try
    files = dir(fullfile(dataPath, subjectFolder, '*xsens*mat'));
    filepath = fullfile(files(1).folder, files(1).name);
    disp(['reading ' filepath '...'])
    xsensTmp = load(filepath);
    xsensData = xsensTmp.xsensData;
    disp('--- done reading Xsens ---')
catch
    warning('No xsens.mat file found, continuing without Xsens data');
end
% xsensData <- xsens center of mass, angles etc

%% load gaze data
plDirs = dir([fullfile(dataPath, subjectFolder, '*pupilcloud_outside*')]);
plDirs = {plDirs.name};

if isempty(plDirs)
    warning('No eyetracker folders found, stopping...');
    return
end

plGazeData = [];
plFixData = [];
plEventData = [];
for f = 1:length(plDirs)
    % import gaze and fixation data
    disp(['reading gaze data from ' fullfile(dataPath, subjectFolder, plDirs{f}) '...'])
    
    pupilFolderPath = fullfile(dataPath, subjectFolder, plDirs{f});
    timeAlignedGazeFile = fullfile(pupilFolderPath, 'time_aligned_gaze.csv');
    
    % Check if time_aligned_gaze.csv exists, if not, create it
    if ~exist(timeAlignedGazeFile, 'file')
        fprintf('⚠️  time_aligned_gaze.csv not found, attempting to create it...\n');
        success = create_time_aligned_gaze(pupilFolderPath);
        if ~success
            warning(['Could not create time-aligned gaze data for folder: ' pupilFolderPath]);
            continue; % Skip this folder
        end
    end
    
    try
        newGaze = importGaze(timeAlignedGazeFile);
        newFix = importFixations(fullfile(pupilFolderPath, 'fixations.csv'));
        newEvent = readtable(fullfile(pupilFolderPath, 'events.csv'));
    catch ex
        warning(['missing files in folder:' newline ex.message]);
        return
    end

    plGazeData = vertcat(plGazeData, newGaze);
    plFixData = vertcat(plFixData, newFix);
    plEventData = vertcat(plEventData, newEvent);
end
disp('--- done reading gaze ---')
% plGazeData, plFixData, plEventData <- pupil labs gaze, fixations, events


%% replace lsl streams with offline recordings


%% eeg
streamsMerged = streams;
% eegIdx = ~cellfun(@isempty,regexp(streamNames,'Counter'));
% eegIdx = contains(streamNames, 'Counter', 'IgnoreCase', true);
eegIdx = cellfun(@(x) contains(x.info.name, {'Counter', 'Sync', 'PROX'}, 'IgnoreCase', true), streamsMerged);

% determine beginning and end of recording
% look for counter gaps in LSL stream and determine b & e for each
% segment to preserve timestamps
gapIdx = find(diff(streams{eegIdx}.time_series) > 1);
gapIdx = [0 gapIdx length(streams{eegIdx}.time_series)];
beg = [];
ending = [];
for g = 1:length(gapIdx)-1
    nextBeg = find(amp.data(cntIdx,:) == streams{eegIdx}.time_series(gapIdx(g)+1));
    if ~isempty(find(amp.data(cntIdx,:) == streams{eegIdx}.time_series(gapIdx(g+1))))
        nextEnding = find(amp.data(cntIdx,:) == streams{eegIdx}.time_series(gapIdx(g+1)));
    % if there is no simple solution, go for this one
    else possibleEnds = find(amp.data(cntIdx,:) ~= 0);
        nextEnding = possibleEnds(end);
    end

    beg = [beg nextBeg];
    ending = [ending nextEnding];
end

% replace counters with eeg data
newEEG = [];
for b = 1:length(beg)
    newEEG = [newEEG amp.data(1:end-1, beg(b):ending(b))];
end
streamsMerged{eegIdx}.time_series = newEEG;
disp('--- replaced EEG data ---');


%% ecg
% farosECGIdx = contains(streamNames, 'faros_ecg', 'IgnoreCase', true);
% farosACCIdx = contains(streamNames, 'faros_acc', 'IgnoreCase', true);
farosECGIdx = cellfun(@(x) contains(x.info.name, 'faros_ecg', 'IgnoreCase', true), streamsMerged);
farosACCIdx = cellfun(@(x) contains(x.info.name, 'faros_acc', 'IgnoreCase', true), streamsMerged);

% resample lsl data to faros sampling rate - only if streams have data
for id = 1:2
    if id == 1
        farosIdx = farosECGIdx;
        streamName = 'ECG';
    elseif id == 2
        farosIdx = farosACCIdx;
        streamName = 'ACC';
    end

    if any(farosIdx)
        lslFaros = streams{farosIdx};
        
        % Check if stream has actual data
        if isempty(lslFaros.time_series) || isempty(lslFaros.time_stamps)
            fprintf('⚠️  Faros %s stream found but empty - skipping resampling\n', streamName);
            continue;
        end
        
        if length(lslFaros.time_stamps) < 2
            fprintf('⚠️  Faros %s stream has insufficient data points (%d) - skipping resampling\n', streamName, length(lslFaros.time_stamps));
            continue;
        end
        
        fprintf('✅ Resampling Faros %s stream (%d samples)\n', streamName, size(lslFaros.time_series, 2));
        
        newLSLTime = lslFaros.time_stamps(1):ecgSampRate:lslFaros.time_stamps(end);
        newLSLData = [];
        for ch = 1:size(lslFaros.time_series, 1)
            newChData = interp1(lslFaros.time_stamps, double(lslFaros.time_series(ch,:)), newLSLTime, 'pchip');
            newLSLData = [newLSLData; newChData];
        end
        streamsMerged{farosIdx}.time_stamps = newLSLTime;
        streamsMerged{farosIdx}.time_series = newLSLData;
    else
        fprintf('ℹ️  No Faros %s stream found - continuing without ECG data\n', streamName);
    end
end
% Process Faros data files only if ECG stream has data
fDataMerged = [];
if any(farosECGIdx) && ~isempty(streamsMerged{farosECGIdx}.time_series) && ~isempty(farosData)
    fprintf('🔗 Processing %d Faros data files for alignment\n', length(farosData));
    
    for f = 1:length(farosData)
        fData = farosData{f};

        % align data from faros and lsl through cross correlation
        try
            [r,lags] = xcorr(streamsMerged{farosECGIdx}.time_series, fData(:,2));
            [~, maxIdx] = max(r(lags >= 0));
            maxIdx = maxIdx + find(lags == 0);
            offset = lags(maxIdx);
            offsetTime = offset*ecgSampRate;

            fData(:,1) = fData(:,1) + streamsMerged{farosECGIdx}.time_stamps(1) + offsetTime;

            fDataMerged = [fDataMerged; fData];
        catch ME
            fprintf('⚠️  Cross-correlation failed for Faros file %d: %s\n', f, ME.message);
            % Continue with raw timestamps if cross-correlation fails
            fDataMerged = [fDataMerged; fData];
        end
    end

    % replace ecg - only if we have valid ECG stream
    if any(farosECGIdx) && ~isempty(streamsMerged{farosECGIdx}.time_stamps) && ~isempty(fDataMerged)
        idx = find(fDataMerged(:,1) >= streamsMerged{farosECGIdx}.time_stamps(1) & fDataMerged(:,1) <= streamsMerged{farosECGIdx}.time_stamps(end));
        if ~isempty(idx)
            streamsMerged{farosECGIdx}.time_stamps = fDataMerged(idx,1)';
            streamsMerged{farosECGIdx}.time_series = fDataMerged(idx,2)';
            fprintf('✅ Replaced ECG data with %d aligned samples\n', length(idx));
        else
            fprintf('⚠️  No overlapping time range found for ECG replacement\n');
        end
    end

    % replace acc - only if we have valid ACC stream
    if any(farosACCIdx) && ~isempty(streamsMerged{farosACCIdx}.time_stamps) && ~isempty(fDataMerged)
        idx = find(fDataMerged(:,1) >= streamsMerged{farosACCIdx}.time_stamps(1) & fDataMerged(:,1) <= streamsMerged{farosACCIdx}.time_stamps(end));
        if ~isempty(idx) && size(fDataMerged, 2) >= 5
            streamsMerged{farosACCIdx}.time_stamps = fDataMerged(idx,1)';
            streamsMerged{farosACCIdx}.time_series = fDataMerged(idx,3:5)';
            fprintf('✅ Replaced ACC data with %d aligned samples\n', length(idx));
        else
            fprintf('⚠️  No overlapping time range or insufficient columns for ACC replacement\n');
        end
    end
    
    disp('--- ECG processing complete ---');
else
    fprintf('ℹ️  Skipping Faros data processing (no ECG stream or no Faros files)\n');
end


%% gaze
% plGazeIdx = contains(streamNames, 'pupil_labs_Gaze', 'IgnoreCase', true);
% plEventIdx = contains(streamNames, 'pupil_labs_Event', 'IgnoreCase', true);
plGazeIdx = cellfun(@(x) contains(x.info.name, 'pupil_labs_Gaze', 'IgnoreCase', true), streamsMerged);
plEventIdx = cellfun(@(x) contains(x.info.name, 'pupil_labs_Event', 'IgnoreCase', true), streamsMerged);

% Check if we have Pupil Labs data available
hasPupilGaze = any(plGazeIdx) && ~isempty(streamsMerged{plGazeIdx}.time_stamps);
hasPupilEvents = any(plEventIdx) && ~isempty(streamsMerged{plEventIdx}.time_stamps);
hasPupilCSV = exist('plGazeData', 'var') && exist('plFixData', 'var') && exist('plEventData', 'var');

if hasPupilCSV && hasPupilGaze
    fprintf('🔗 Processing Pupil Labs data integration\n');
    
    % add lsl times to fixations
    if size(plFixData, 1) > 0
        plFixData.lsl_times = zeros(size(plFixData,1),1);

        for i = 1:size(plFixData,1)
            fixNum = plFixData.fixationId(i);
            fixId = find(plGazeData.fixationId == fixNum);
            if ~isempty(fixId)
                fixId = fixId(1); % use start of fixation
                plFixData.lsl_times(i) = plGazeData.lsl_times(fixId);
            end
        end
        fprintf('✅ Added LSL times to %d fixations\n', size(plFixData,1));
    else
        fprintf('ℹ️  No fixation data available\n');
    end

    % add lsl times to events
    if size(plEventData, 1) > 0 && hasPupilEvents
        plEventData.lsl_times = zeros(size(plEventData,1),1);

        % add known event times from LSL data
        for e = 1:length(streamsMerged{plEventIdx}.time_stamps)
            idx = find(strcmp(plEventData.name, streamsMerged{plEventIdx}.time_series(e)));
            if ~isempty(idx)
                plEventData.lsl_times(idx) = streamsMerged{plEventIdx}.time_stamps(e);
            end
        end

        % calculate times for remaining events from nearest sync point
        for e = 1:size(plEventData,1)
            if plEventData.lsl_times(e) == 0
                idx = find(plEventData.lsl_times ~= 0);
                if ~isempty(idx)
                    [~,refIdx] = min(abs(idx - e));
                    plEventData.lsl_times(e) = plEventData.lsl_times(idx(refIdx)) + (plEventData.timestamp_ns_(e) - plEventData.timestamp_ns_(idx(refIdx)))/10^9;
                end
            end
        end
        fprintf('✅ Added LSL times to %d events\n', size(plEventData,1));
    else
        fprintf('ℹ️  No event data available or no event stream\n');
    end

    % replace lsl data with pupil cloud data if we have valid gaze data
    if size(plGazeData, 1) > 0
        % gaze:
        gazeIdx = find(plGazeData.lsl_times >= streamsMerged{plGazeIdx}.time_stamps(1) & plGazeData.lsl_times <= streamsMerged{plGazeIdx}.time_stamps(end));
        
        if ~isempty(gazeIdx)
            gazeDataNew = plGazeData(gazeIdx,:);
            streamsMerged{plGazeIdx}.time_stamps = gazeDataNew.lsl_times';
            streamsMerged{plGazeIdx}.time_series = [gazeDataNew.gazeXpx, gazeDataNew.gazeYpx, gazeDataNew.fixationId, gazeDataNew.blinkId, gazeDataNew.azimuthdeg, gazeDataNew.elevationdeg]';

            streamsMerged{plGazeIdx}.info.desc.channels.channel{3}.label = 'fixationId';
            streamsMerged{plGazeIdx}.info.desc.channels.channel{3}.eye = 'both';
            streamsMerged{plGazeIdx}.info.desc.channels.channel{4}.label = 'blinkId';
            streamsMerged{plGazeIdx}.info.desc.channels.channel{4}.eye = 'both';
            streamsMerged{plGazeIdx}.info.desc.channels.channel{5}.label = 'azimuthdeg';
            streamsMerged{plGazeIdx}.info.desc.channels.channel{5}.eye = 'both';
            streamsMerged{plGazeIdx}.info.desc.channels.channel{6}.label = 'elevationdeg';
            streamsMerged{plGazeIdx}.info.desc.channels.channel{6}.eye = 'both';
            
            fprintf('✅ Replaced gaze data with %d samples from CSV\n', length(gazeIdx));
        else
            fprintf('⚠️  No overlapping gaze data found\n');
        end

        % fixations (new channel):
        if size(plFixData, 1) > 0 && any(isfield(plFixData, 'lsl_times')) && any(plFixData.lsl_times > 0)
            fixIdx = find(plFixData.lsl_times >= streamsMerged{plGazeIdx}.time_stamps(1) & plFixData.lsl_times <= streamsMerged{plGazeIdx}.time_stamps(end));
            
            if ~isempty(fixIdx)
                fixDataNew = plFixData(fixIdx,:);
                plFixIdx = length(streamsMerged) + 1;
                streamsMerged{plFixIdx}.time_stamps = fixDataNew.lsl_times';
                streamsMerged{plFixIdx}.time_series = [fixDataNew.fixationId, fixDataNew.durationms, fixDataNew.fixationXpx, fixDataNew.fixationYpx, fixDataNew.azimuthdeg, fixDataNew.elevationdeg]';
                fprintf('✅ Added fixation stream with %d samples\n', length(fixIdx));
            else
                fprintf('ℹ️  No overlapping fixation data found\n');
            end
        else
            fprintf('ℹ️  No valid fixation data to add\n');
        end
    else
        fprintf('ℹ️  No gaze data available\n');
    end
    
        % Add fixation stream metadata if it was created
        if exist('plFixIdx', 'var')
            streamsMerged{plFixIdx}.info.name = 'pupil_labs_Fixations';
            streamsMerged{plFixIdx}.info.type = 'fixations';
            streamsMerged{plFixIdx}.info.channel_count = '6';
            streamsMerged{plFixIdx}.info.desc.channels.channel{1}.label = 'fixationId';
            streamsMerged{plFixIdx}.info.desc.channels.channel{1}.eye = 'both';
            streamsMerged{plFixIdx}.info.desc.channels.channel{2}.label = 'durationMs';
            streamsMerged{plFixIdx}.info.desc.channels.channel{2}.eye = 'both';
            streamsMerged{plFixIdx}.info.desc.channels.channel{3}.label = 'fixationX';
            streamsMerged{plFixIdx}.info.desc.channels.channel{3}.eye = 'both';
            streamsMerged{plFixIdx}.info.desc.channels.channel{4}.label = 'fixationY';
            streamsMerged{plFixIdx}.info.desc.channels.channel{4}.eye = 'both';
            streamsMerged{plFixIdx}.info.desc.channels.channel{5}.label = 'azimuthDeg';
            streamsMerged{plFixIdx}.info.desc.channels.channel{5}.eye = 'both';
            streamsMerged{plFixIdx}.info.desc.channels.channel{6}.label = 'elevationDeg';
            streamsMerged{plFixIdx}.info.desc.channels.channel{6}.eye = 'both';
        end

        % Update events if we have event data
        if hasPupilEvents && size(plEventData, 1) > 0
            streamsMerged{plEventIdx}.time_stamps = plEventData.lsl_times';
            streamsMerged{plEventIdx}.time_series = horzcat(plEventData.name, plEventData.type)';
            streamsMerged{plEventIdx}.info.channel_count = '2';
            fprintf('✅ Updated event stream with %d events\n', size(plEventData, 1));
        end
    
    disp('--- Pupil Labs processing complete ---');
else
    if ~hasPupilCSV
        fprintf('ℹ️  No Pupil Labs CSV data found - skipping Pupil integration\n');
    elseif ~hasPupilGaze
        fprintf('ℹ️  No Pupil Labs LSL stream found - skipping Pupil integration\n');
    end
end


%% xsens
xsensCoMIdx = cellfun(@(x) contains(x.info.name, 'CenterOfMass1', 'IgnoreCase', true), streamsMerged);
xsensSegmIdx = cellfun(@(x) contains(x.info.name, 'LinearSegmentKinematicsDatagram1', 'IgnoreCase', true), streamsMerged);
xsensJointsIdx = cellfun(@(x) contains(x.info.name, 'JointAnglesDatagram1', 'IgnoreCase', true), streamsMerged);
xsensEulerIdx = cellfun(@(x) contains(x.info.name, 'EulerDatagram1', 'IgnoreCase', true), streamsMerged);
xsensQuatIdx = cellfun(@(x) contains(x.info.name, 'QuaternionDatagram1', 'IgnoreCase', true), streamsMerged);

% Check XSens data and stream availability
hasXSensData = ~isempty(xsensData);
hasXSensCoMStream = any(xsensCoMIdx) && ~isempty(streams{xsensCoMIdx}.time_stamps);

if hasXSensData && hasXSensCoMStream
    fprintf('🔗 Processing XSens data integration\n');
    
    xsensLSLTime = [];
    xsensCoMData = [];
    xsensSegmentData = [];
    xsensJointData = [];
    xsensEulerData = [];
    xsensQuatData = [];
    
    try
        % define vars
        frameRate = xsensData.info{1}.MVN_version(2);
        comData = xsensData.center_of_mass{1};
        
        % Check if streams have sufficient data for cross-correlation
        if length(streams{xsensCoMIdx}.time_stamps) < 2
            fprintf('⚠️  Insufficient XSens LSL data for alignment - using raw timestamps\n');
            xsensTime = (0:(size(comData,1)-1)) * (1/frameRate);
        else
            xsensTime = (0:(size(comData,1)-1)) * (1/frameRate) + streams{xsensCoMIdx}.time_stamps(1);

            % resample existing LSL data to remove time gaps
            lslCoMX = streams{xsensCoMIdx}.time_series(1,:);
            newLSLTime = streams{xsensCoMIdx}.time_stamps(1):1/frameRate:streams{xsensCoMIdx}.time_stamps(end);
            
            if length(newLSLTime) > 1
                newCoMX = interp1(streams{xsensCoMIdx}.time_stamps, lslCoMX, newLSLTime, 'linear', 'extrap');

                % use correlation to find offset between offline and LSL data
                if length(newCoMX) > 1 && length(comData.CoM_pos_x) > 1
                    [r,lags] = xcorr(newCoMX, comData.CoM_pos_x);
                    [~,maxIdx] = max(r);
                    xsensTime = xsensTime + lags(maxIdx)*(1/frameRate);
                    fprintf('✅ XSens alignment completed with cross-correlation\n');
                else
                    fprintf('⚠️  Insufficient data for cross-correlation - using raw timestamps\n');
                end
            else
                fprintf('⚠️  Insufficient resampled data - using raw timestamps\n');
            end
        end

        % identify new time array
        if exist('newLSLTime', 'var') && ~isempty(newLSLTime)
            idx = find(xsensTime >= newLSLTime(1) & xsensTime <= newLSLTime(end));
        else
            % Use all available data if no LSL time reference
            idx = 1:length(xsensTime);
        end
        
        if ~isempty(idx)
            xsensLSLTime = [xsensLSLTime xsensTime(idx)];

            % copy center of mass data
            if isfield(xsensData, 'center_of_mass') && ~isempty(xsensData.center_of_mass{1})
                for x = 2:size(xsensData.center_of_mass{1},2)
                    tmp = table2array(xsensData.center_of_mass{1});
                    tmp = tmp(:,x);
                    xsensCoMData = [xsensCoMData; tmp(idx,:)'];
                end
                fprintf('✅ Processed %d center of mass channels\n', size(xsensData.center_of_mass{1},2)-1);
            end

            % copy segment positions
            if isfield(xsensData, 'segment_position') && ~isempty(xsensData.segment_position{1})
                for x = 2:size(xsensData.segment_position{1},2)
                    tmp = table2array(xsensData.segment_position{1});
                    tmp = tmp(:,x);
                    xsensSegmentData = [xsensSegmentData; tmp(idx,:)'];
                end
                fprintf('✅ Processed %d segment position channels\n', size(xsensData.segment_position{1},2)-1);
            end

            % copy joint angles
            if isfield(xsensData, 'joint_angles') && ~isempty(xsensData.joint_angles{1})
                for x = 2:size(xsensData.joint_angles{1},2)
                    tmp = table2array(xsensData.joint_angles{1});
                    tmp = tmp(:,x);
                    xsensJointData = [xsensJointData; tmp(idx,:)'];
                end
                fprintf('✅ Processed %d joint angle channels\n', size(xsensData.joint_angles{1},2)-1);
            end

            % copy euler orientation
            if isfield(xsensData, 'segment_orientation_euler') && ~isempty(xsensData.segment_orientation_euler{1})
                for x = 2:size(xsensData.segment_orientation_euler{1},2)
                    tmp = table2array(xsensData.segment_orientation_euler{1});
                    tmp = tmp(:,x);
                    xsensEulerData = [xsensEulerData; tmp(idx,:)'];
                end
                fprintf('✅ Processed %d Euler orientation channels\n', size(xsensData.segment_orientation_euler{1},2)-1);
            end

            % copy quaternions
            if isfield(xsensData, 'segment_orientation_quat') && ~isempty(xsensData.segment_orientation_quat{1})
                for x = 1:size(xsensData.segment_orientation_quat{1},2)
                    tmp = table2array(xsensData.segment_orientation_quat{1});
                    tmp = tmp(:,x);
                    xsensQuatData = [xsensQuatData; tmp(idx,:)'];
                end
                fprintf('✅ Processed %d quaternion channels\n', size(xsensData.segment_orientation_quat{1},2));
            end
        else
            fprintf('⚠️  No overlapping time indices found between XSens and LSL data\n');
        end

        % write data to xdf streams (only if streams exist and data was processed)
        if ~isempty(xsensCoMData) && any(xsensCoMIdx)
            streamsMerged{xsensCoMIdx}.time_stamps = xsensLSLTime;
            streamsMerged{xsensCoMIdx}.time_series = xsensCoMData;
            streamsMerged{xsensCoMIdx}.info.desc.channels = {xsensData.center_of_mass{1}.Properties.VariableNames{2:end}};
        end

        % write segment positions to xdf streams
        if ~isempty(xsensSegmentData) && any(xsensSegmIdx)
            streamsMerged{xsensSegmIdx}.time_stamps = xsensLSLTime;
            streamsMerged{xsensSegmIdx}.time_series = xsensSegmentData;
            streamsMerged{xsensSegmIdx}.info.desc.channels = {xsensData.segment_position{1}.Properties.VariableNames{2:end}};
        end

        % write joint angles to xdf streams
        if ~isempty(xsensJointData) && any(xsensJointsIdx)
            streamsMerged{xsensJointsIdx}.time_stamps = xsensLSLTime;
            streamsMerged{xsensJointsIdx}.time_series = xsensJointData;
            streamsMerged{xsensJointsIdx}.info.desc.channels = {xsensData.joint_angles{1}.Properties.VariableNames{2:end}};
        end

        % write euler to xdf streams
        if ~isempty(xsensEulerData) && any(xsensEulerIdx)
            streamsMerged{xsensEulerIdx}.time_stamps = xsensLSLTime;
            streamsMerged{xsensEulerIdx}.time_series = xsensEulerData;
            streamsMerged{xsensEulerIdx}.info.desc.channels = {xsensData.segment_orientation_euler{1}.Properties.VariableNames{2:end}};
        end

        % write quaternions to xdf streams
        if ~isempty(xsensQuatData) && any(xsensQuatIdx)
            streamsMerged{xsensQuatIdx}.time_stamps = xsensLSLTime;
            streamsMerged{xsensQuatIdx}.time_series = xsensQuatData;
            streamsMerged{xsensQuatIdx}.info.desc.channels = {xsensData.segment_orientation_quat{1}.Properties.VariableNames{2:end}};
        end
        
        disp('--- XSens processing complete ---');
        
    catch ME
        fprintf('⚠️  XSens data processing failed: %s\n', ME.message);
        disp('--- XSens processing failed, continuing without XSens data ---');
    end
else
    if ~hasXSensData
        fprintf('ℹ️  No XSens data available - skipping XSens integration\n');
    elseif ~hasXSensCoMStream
        fprintf('ℹ️  No XSens LSL stream available - skipping XSens integration\n');
    end
end


%% add lsl data to eeg struct

fprintf('🔗 Integrating LSL streams into EEG structure\n');

% match counter values in eeg and lsl data
[~, idxAmp, idxLSL] = intersect(amp.data(end,:), streams{eegIdx}.time_series);
counterLSLTime = streams{eegIdx}.time_stamps(idxLSL);

fprintf('✅ Found %d matching counter values for synchronization\n', length(idxAmp));

newChannelIdx = amp.nbchan+1;

% sample all lsl data at the times of the counter values before adding

% add ecg
if any(farosECGIdx) && ~isempty(streamsMerged{farosECGIdx}.time_stamps) && ~isempty(streamsMerged{farosECGIdx}.time_series)
    try
        ampECG = interp1(streamsMerged{farosECGIdx}.time_stamps, streamsMerged{farosECGIdx}.time_series, counterLSLTime, 'cubic');
        amp.data(newChannelIdx,:) = NaN(1, size(amp.data, 2));
        amp.data(newChannelIdx,idxAmp) = ampECG;
        amp.chanlocs(newChannelIdx).labels = 'Faros ECG';
        amp.chanlocs(newChannelIdx).type = 'ECG';
        newChannelIdx = newChannelIdx + 1;
        fprintf('✅ Merged ECG into EEG structure\n');
    catch ME
        fprintf('⚠️  Failed to merge ECG data: %s\n', ME.message);
    end
else
    fprintf('ℹ️  No ECG data available for merging\n');
end

% add photo
% photoIdx = contains(streamNames, 'Photo', 'IgnoreCase', true);
photoIdx = cellfun(@(x) contains(x.info.name, 'Photo', 'IgnoreCase', true), streamsMerged);
if any(photoIdx) && ~isempty(streamsMerged{photoIdx}.time_stamps) && ~isempty(streamsMerged{photoIdx}.time_series)
    try
        ampPhoto = interp1(streamsMerged{photoIdx}.time_stamps, double(streamsMerged{photoIdx}.time_series), counterLSLTime, 'pchip');
        amp.data(newChannelIdx,:) = NaN(1, size(amp.data, 2));
        amp.data(newChannelIdx, idxAmp) = ampPhoto;
        amp.chanlocs(newChannelIdx).labels = 'PhotoSensor';
        amp.chanlocs(newChannelIdx).type = 'Photo';
        newChannelIdx = newChannelIdx + 1;
        fprintf('✅ Merged PhotoSensor into EEG structure\n');
    catch ME
        fprintf('⚠️  Failed to merge PhotoSensor data: %s\n', ME.message);
    end
else
    fprintf('ℹ️  No PhotoSensor data available for merging\n');
end

% add gaze
if any(plGazeIdx) && ~isempty(streamsMerged{plGazeIdx}.time_stamps) && ~isempty(streamsMerged{plGazeIdx}.time_series)
    try
        % Check that we have enough channels in the gaze data
        if size(streamsMerged{plGazeIdx}.time_series, 1) >= 4
            ampGazeX = interp1(streamsMerged{plGazeIdx}.time_stamps, streamsMerged{plGazeIdx}.time_series(1,:), counterLSLTime, 'pchip');
            ampGazeY = interp1(streamsMerged{plGazeIdx}.time_stamps, streamsMerged{plGazeIdx}.time_series(2,:), counterLSLTime, 'pchip');
            ampFixID = interp1(streamsMerged{plGazeIdx}.time_stamps, streamsMerged{plGazeIdx}.time_series(3,:), counterLSLTime, 'nearest');
            ampBlinkID = interp1(streamsMerged{plGazeIdx}.time_stamps, streamsMerged{plGazeIdx}.time_series(4,:), counterLSLTime, 'nearest');
            amp.data(newChannelIdx:newChannelIdx+3,:) = NaN(4, size(amp.data, 2));
            amp.data(newChannelIdx:newChannelIdx+3,idxAmp) = [ampGazeX; ampGazeY; ampFixID; ampBlinkID];
            for c = 0:3
                amp.chanlocs(newChannelIdx+c).type = 'Gaze';
            end
            amp.chanlocs(newChannelIdx).labels = 'GazeX';
            amp.chanlocs(newChannelIdx+1).labels = 'GazeY';
            amp.chanlocs(newChannelIdx+2).labels = 'FixationID';
            amp.chanlocs(newChannelIdx+3).labels = 'BlinkID';
            newChannelIdx = newChannelIdx + 4;
            fprintf('✅ Merged gaze data into EEG structure (4 channels)\n');
        else
            fprintf('⚠️  Insufficient gaze channels (%d < 4) - skipping gaze integration\n', size(streamsMerged{plGazeIdx}.time_series, 1));
        end
    catch ME
        fprintf('⚠️  Failed to merge gaze data: %s\n', ME.message);
    end
else
    fprintf('ℹ️  No gaze data available for merging\n');
end

%% add xsens
if ~isempty(xsensData) && any(xsensCoMIdx) && ~isempty(streamsMerged{xsensCoMIdx}.time_stamps)
    try
        fprintf('🔗 Integrating XSens data into EEG structure\n');
        
        % Center of mass (keep existing)
        for m = 1:size(streamsMerged{xsensCoMIdx}.time_series,1)
            ampCoM = interp1(streamsMerged{xsensCoMIdx}.time_stamps, streamsMerged{xsensCoMIdx}.time_series(m,:), counterLSLTime, 'pchip');
            amp.data(newChannelIdx+m-1,:) = NaN(1, size(amp.data, 2));
            amp.data(newChannelIdx+m-1, idxAmp) = ampCoM; % Fixed: was ampPhoto
            amp.chanlocs(newChannelIdx+m-1).labels = char(streamsMerged{xsensCoMIdx}.info.desc.channels(m));
            amp.chanlocs(newChannelIdx+m-1).type = 'Xsens';
        end
        newChannelIdx = newChannelIdx + size(streamsMerged{xsensCoMIdx}.time_series,1);
        
        fprintf('✅ Added %d XSens CoM channels to EEG\n', size(streamsMerged{xsensCoMIdx}.time_series,1));
        
        % Critical gait joints - only if joint data is available
        if any(xsensJointsIdx) && ~isempty(streamsMerged{xsensJointsIdx}.time_stamps)
            % Ankle dorsiflexion/plantarflexion (MOST IMPORTANT)
            ankleChannels = {'Right_Ankle_Dorsiflexion_Plantarflexion', 'Left_Ankle_Dorsiflexion_Plantarflexion'};
            for ac = 1:length(ankleChannels)
                ankleIdx = contains(streamsMerged{xsensJointsIdx}.info.desc.channels, ankleChannels{ac}, 'IgnoreCase', true);
                if any(ankleIdx)
                    ampAnkleData = interp1(streamsMerged{xsensJointsIdx}.time_stamps, streamsMerged{xsensJointsIdx}.time_series(ankleIdx,:), counterLSLTime, 'pchip');
                    amp.data(newChannelIdx,:) = NaN(1, size(amp.data, 2));
                    amp.data(newChannelIdx, idxAmp) = ampAnkleData;
                    amp.chanlocs(newChannelIdx).labels = ankleChannels{ac};
                    amp.chanlocs(newChannelIdx).type = 'Xsens_Ankle';
                    newChannelIdx = newChannelIdx + 1;
                end
            end
            
            % Hip flexion/extension (important for swing/stance detection)
            hipChannels = {'Right_Hip_Flexion_Extension', 'Left_Hip_Flexion_Extension'};
            for hc = 1:length(hipChannels)
                hipIdx = contains(streamsMerged{xsensJointsIdx}.info.desc.channels, hipChannels{hc}, 'IgnoreCase', true);
                if any(hipIdx)
                    ampHipData = interp1(streamsMerged{xsensJointsIdx}.time_stamps, streamsMerged{xsensJointsIdx}.time_series(hipIdx,:), counterLSLTime, 'pchip');
                    amp.data(newChannelIdx,:) = NaN(1, size(amp.data, 2));
                    amp.data(newChannelIdx, idxAmp) = ampHipData;
                    amp.chanlocs(newChannelIdx).labels = hipChannels{hc};
                    amp.chanlocs(newChannelIdx).type = 'Xsens_Hip';
                    newChannelIdx = newChannelIdx + 1;
                end
            end
            
            % Knee flexion/extension (important for gait cycle)
            kneeChannels = {'Right_Knee_Flexion_Extension', 'Left_Knee_Flexion_Extension'};
            for kc = 1:length(kneeChannels)
                kneeIdx = contains(streamsMerged{xsensJointsIdx}.info.desc.channels, kneeChannels{kc}, 'IgnoreCase', true);
                if any(kneeIdx)
                    ampKneeData = interp1(streamsMerged{xsensJointsIdx}.time_stamps, streamsMerged{xsensJointsIdx}.time_series(kneeIdx,:), counterLSLTime, 'pchip');
                    amp.data(newChannelIdx,:) = NaN(1, size(amp.data, 2));
                    amp.data(newChannelIdx, idxAmp) = ampKneeData;
                    amp.chanlocs(newChannelIdx).labels = kneeChannels{kc};
                    amp.chanlocs(newChannelIdx).type = 'Xsens_Knee';
                    newChannelIdx = newChannelIdx + 1;
                end
            end
            
            % Ball of foot (toe) flexion/extension (important for toe-off detection)
            toeChannels = {'Right_Ball_Foot_Flexion_Extension', 'Left_Ball_Foot_Flexion_Extension'};
            for tc = 1:length(toeChannels)
                toeIdx = contains(streamsMerged{xsensJointsIdx}.info.desc.channels, toeChannels{tc}, 'IgnoreCase', true);
                if any(toeIdx)
                    ampToeData = interp1(streamsMerged{xsensJointsIdx}.time_stamps, streamsMerged{xsensJointsIdx}.time_series(toeIdx,:), counterLSLTime, 'pchip');
                    amp.data(newChannelIdx,:) = NaN(1, size(amp.data, 2));
                    amp.data(newChannelIdx, idxAmp) = ampToeData;
                    amp.chanlocs(newChannelIdx).labels = toeChannels{tc};
                    amp.chanlocs(newChannelIdx).type = 'Xsens_Toe';
                    newChannelIdx = newChannelIdx + 1;
                end
            end
            fprintf('✅ Added XSens joint angle channels for gait analysis\n');
        end
        
        % Foot positions - only if segment data is available
        if any(xsensSegmIdx) && ~isempty(streamsMerged{xsensSegmIdx}.time_stamps)
            footChannels = {'Left_Foot_Z', 'Right_Foot_Z'}; % Z for vertical movement
            for fc = 1:length(footChannels)
                footIdx = contains(streamsMerged{xsensSegmIdx}.info.desc.channels, footChannels{fc}, 'IgnoreCase', true);
                if any(footIdx)
                    ampFootData = interp1(streamsMerged{xsensSegmIdx}.time_stamps, streamsMerged{xsensSegmIdx}.time_series(footIdx,:), counterLSLTime, 'pchip');
                    amp.data(newChannelIdx,:) = NaN(1, size(amp.data, 2));
                    amp.data(newChannelIdx, idxAmp) = ampFootData;
                    amp.chanlocs(newChannelIdx).labels = footChannels{fc};
                    amp.chanlocs(newChannelIdx).type = 'Xsens_Foot_Pos';
                    newChannelIdx = newChannelIdx + 1;
                end
            end
            fprintf('✅ Added XSens foot position channels for gait analysis\n');
        end
        
        disp('--- XSens integration into EEG complete ---');
        
    catch ME
        fprintf('⚠️  XSens EEG integration failed: %s\n', ME.message);
    end
else
    fprintf('ℹ️  No XSens data available for EEG integration\n');
end


%% Enhanced Pupil Event Processing from CSV files
% Load and integrate pupil events directly from CSV files with robust time synchronization
fprintf('🔄 Loading enhanced pupil events from CSV files...\n');

% Initialize pupil data structure
pupilData = struct();
pupilData.events = table();
pupilData.gaze = table();
pupilData.fixations = table();
pupilData.blinks = table();

% Process all pupil folders with comprehensive CSV loading
for f = 1:length(plDirs)
    pupilDataPath = fullfile(dataPath, subjectFolder, plDirs{f});
    
    fprintf('   📁 Processing folder: %s\n', plDirs{f});
    
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
                    fprintf('     ⚠️ %s too small (%d bytes), skipping\n', csvFiles{i}, fileInfo.bytes);
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
                    fprintf('     📥 Loaded %s: %d rows\n', csvFiles{i}, height(data));
                else
                    fprintf('     ⚠️ %s is empty\n', csvFiles{i});
                end
            catch ME
                fprintf('     ❌ Could not load %s: %s\n', csvFiles{i}, ME.message);
            end
        else
            if i <= 2 % Only warn for critical files
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
                    fprintf('     🔗 Merged %s data (%d new rows)\n', dataTypeStr, height(currentPupilData.(dataTypeStr)));
                catch ME
                    fprintf('     ⚠️ Could not merge %s data: %s\n', dataTypeStr, ME.message);
                end
            end
        end
    end
end

% Check if we have required data
if ~isfield(pupilData, 'events') | height(pupilData.events) == 0
    fprintf('   ⚠️ No pupil events found for synchronization\n');
    return;
end

% Summary of loaded Pupil Labs data
fprintf('   📊 PUPIL LABS DATA SUMMARY:\n');
dataTypes = {'events', 'gaze', 'fixations', 'blinks'};
for dataType = dataTypes
    if isfield(pupilData, dataType{1}) & height(pupilData.(dataType{1})) > 0
        fprintf('     - %s: %d rows\n', dataType{1}, height(pupilData.(dataType{1})));
    else
        fprintf('     - %s: not available\n', dataType{1});
    end
end

%% Robust Time Synchronization
fprintf('   🔄 Performing robust time alignment...\n');

% Extract sync events from XDF (try multiple methods)
xdfSyncEvents = [];
xdfSyncTimes = [];

% Find marker/event streams
for i = 1:length(streamsMerged)
    if isfield(streamsMerged{i}, 'time_series') & isfield(streamsMerged{i}, 'info') && ...
       (isfield(streamsMerged{i}.info, 'type') & ...
        (strcmpi(streamsMerged{i}.info.type, 'Markers') | strcmpi(streamsMerged{i}.info.type, 'Event')))
        
        events = streamsMerged{i}.time_series(1,:);
        times = streamsMerged{i}.time_stamps;
        
        if iscell(events)
            % Method 1: Look for lsl.time_sync events (preferred)
            syncMask = contains(events, 'lsl.time_sync');
            if any(syncMask)
                xdfSyncEvents = [xdfSyncEvents; events(syncMask)];
                xdfSyncTimes = [xdfSyncTimes; times(syncMask)'];
                fprintf('      ✅ Found %d lsl.time_sync events\n', sum(syncMask));
            end
            
            % Method 2: Look for recording events (fallback)
            recordingMask = contains(events, 'recording.begin') | contains(events, 'recording.end');
            if any(recordingMask) & length(xdfSyncEvents) < 2
                xdfSyncEvents = [xdfSyncEvents; events(recordingMask)];
                xdfSyncTimes = [xdfSyncTimes; times(recordingMask)'];
                fprintf('      ✅ Found %d recording.begin/end events\n', sum(recordingMask));
            end
        end
    end
end

% Determine sync method
lslSyncCount = sum(contains(xdfSyncEvents, 'lsl.time_sync'));
recordingCount = sum(contains(xdfSyncEvents, 'recording.'));

if lslSyncCount >= 2
    syncMethod = 'lsl.time_sync';
    fprintf('   🎯 Using LSL time sync events (%d found)\n', lslSyncCount);
elseif recordingCount >= 2
    syncMethod = 'recording.begin/end';
    fprintf('   🎯 Using recording begin/end events (%d found)\n', recordingCount);
else
    syncMethod = 'insufficient';
    fprintf('   ❌ Insufficient sync events found (LSL: %d, Recording: %d)\n', lslSyncCount, recordingCount);
    return;
end

% Find timestamp column in pupil CSV
timestampCol = find_timestamp_column(pupilData.events);

% Extract corresponding sync events from CSV
csvSyncEvents = [];
csvSyncTimesNs = [];

% Match sync events
for i = 1:length(xdfSyncEvents)
    eventName = xdfSyncEvents{i};
    csvMatch = strcmp(pupilData.events.name, eventName);
    if any(csvMatch)
        csvSyncEvents = [csvSyncEvents; {eventName}];
        matchIdx = find(csvMatch, 1, 'first');
        csvSyncTimesNs = [csvSyncTimesNs; pupilData.events.(timestampCol)(matchIdx)];
    end
end

if length(csvSyncEvents) < 2
    fprintf('   ❌ Insufficient matching sync events (XDF: %d, CSV: %d)\n', length(xdfSyncEvents), length(csvSyncEvents));
    return;
end

% Sort events by time
[xdfSyncTimes, xdfSortIdx] = sort(xdfSyncTimes);
xdfSyncEvents = xdfSyncEvents(xdfSortIdx);

[csvSyncTimesNs, csvSortIdx] = sort(csvSyncTimesNs);
csvSyncEvents = csvSyncEvents(csvSortIdx);

% Convert CSV times to seconds
csvSyncSeconds = csvSyncTimesNs * 1e-9;

% ROBUST TIME ALIGNMENT using time differences to avoid numerical issues
if length(xdfSyncTimes) >= 2
    % Method: Use time differences between consecutive sync events
    xdfDiffs = diff(xdfSyncTimes);
    csvDiffs = diff(csvSyncSeconds);
    
    fprintf('   🔧 Using time differences method for alignment\n');
    fprintf('      XDF time diffs: %.3f ± %.3f seconds\n', mean(xdfDiffs), std(xdfDiffs));
    fprintf('      CSV time diffs: %.3f ± %.3f seconds\n', mean(csvDiffs), std(csvDiffs));
    
    % Calculate slope (should be close to 1.0 for good sync)
    if length(xdfDiffs) > 1
        slope = xdfDiffs \ csvDiffs;  % Robust slope estimation
    else
        slope = xdfDiffs / csvDiffs;
    end
    
    % Calculate intercept using the first matched pair
    intercept = xdfSyncTimes(1) - slope * csvSyncSeconds(1);
    
    % Quality metrics
    predicted_xdf = intercept + slope * csvSyncSeconds;
    residuals = xdfSyncTimes - predicted_xdf;
    rmse = sqrt(mean(residuals.^2));
    r_squared = 1 - sum(residuals.^2) / sum((xdfSyncTimes - mean(xdfSyncTimes)).^2);
    
    fprintf('   ✅ Time alignment: intercept=%.6f, slope=%.9f\n', intercept, slope);
    fprintf('   📊 Quality: R²=%.6f, RMSE=%.6f sec (%d events)\n', r_squared, rmse, length(xdfSyncTimes));
    
    % Validate alignment quality
    if abs(slope - 1.0) > 0.1
        warning('Clock drift detected: slope = %.6f', slope);
    end
    
    if rmse > 0.5
        warning('High RMSE (%.3f) - alignment may be poor', rmse);
    end
    
    % Convert all pupil event times to XDF time
    csv_to_xdf_time = @(csv_time_ns) intercept + slope * (csv_time_ns * 1e-9);
    pupilData.events.timestamp_xdf = csv_to_xdf_time(pupilData.events.(timestampCol));
    
    % Add comprehensive pupil events to EEG structure
    lslTimeInEEG = amp.times / 1000; % Convert to seconds
    refIdx = find(amp.data(cntIdx,:) == streams{eegIdx}.time_series(1));
    lslTimeInEEG = lslTimeInEEG - lslTimeInEEG(refIdx) + streams{eegIdx}.time_stamps(1);
    
    eventsAdded = 0;
    blinksAdded = 0;
    fixationsAdded = 0;
    currentEventCount = length(amp.event);
    
    % Add events from events.csv
    for i = 1:height(pupilData.events)
        [timeDiff, sampleIdx] = min(abs(lslTimeInEEG - pupilData.events.timestamp_xdf(i)));
        
        if timeDiff < 0.1  % Within 100ms
            currentEventCount = currentEventCount + 1;
            amp.event(currentEventCount).latency = sampleIdx;
            amp.event(currentEventCount).duration = 1;
            amp.event(currentEventCount).type = char(pupilData.events.name(i));
            amp.event(currentEventCount).source = 'pupil_csv_enhanced';
            
            % Add comprehensive metadata
            eventCols = pupilData.events.Properties.VariableNames;
            if ismember('type', eventCols)
                amp.event(currentEventCount).pupil_type = char(pupilData.events.type(i));
            end
            
            amp.event(currentEventCount).pupil_timestamp_ns = pupilData.events.(timestampCol)(i);
            amp.event(currentEventCount).pupil_timestamp_xdf = pupilData.events.timestamp_xdf(i);
            
            % Add section and recording IDs
            if ismember('sectionId', eventCols)
                amp.event(currentEventCount).pupil_section_id = char(pupilData.events.sectionId(i));
            end
            if ismember('recordingId', eventCols)
                amp.event(currentEventCount).pupil_recording_id = char(pupilData.events.recordingId(i));
            end
            
            eventsAdded = eventsAdded + 1;
        end
    end
    
    % Add blinks as events
    if isfield(pupilData, 'blinks') && height(pupilData.blinks) > 0
        % Align blink timestamps
        blinkTimestampCol = find_timestamp_column(pupilData.blinks);
        if ~isempty(blinkTimestampCol)
            pupilData.blinks.timestamp_xdf = csv_to_xdf_time(pupilData.blinks.(blinkTimestampCol));
            
            for i = 1:height(pupilData.blinks)
                [timeDiff, sampleIdx] = min(abs(lslTimeInEEG - pupilData.blinks.timestamp_xdf(i)));
                
                if timeDiff < 0.1
                    currentEventCount = currentEventCount + 1;
                    amp.event(currentEventCount).latency = sampleIdx;
                    amp.event(currentEventCount).type = 'blink';
                    amp.event(currentEventCount).source = 'pupil_blinks';
                    amp.event(currentEventCount).duration = 1;
                    
                    % Add blink duration if available
                    durCols = {'duration_ms_', 'durationms', 'duration'};
                    for col = durCols
                        if ismember(col{1}, pupilData.blinks.Properties.VariableNames)
                            duration_ms = pupilData.blinks.(col{1})(i);
                            if ~isnan(duration_ms)
                                amp.event(currentEventCount).duration = round(duration_ms * amp.srate / 1000);
                                amp.event(currentEventCount).blink_duration_ms = duration_ms;
                            end
                            break;
                        end
                    end
                    
                    blinksAdded = blinksAdded + 1;
                end
            end
        end
    end
    
    % Add fixations as events
    if isfield(pupilData, 'fixations') && height(pupilData.fixations) > 0
        % Find start timestamp column for fixations
        startCol = find(contains(pupilData.fixations.Properties.VariableNames, 'startTimestamp', 'IgnoreCase', true));
        if ~isempty(startCol)
            colName = pupilData.fixations.Properties.VariableNames{startCol(1)};
            pupilData.fixations.timestamp_xdf = csv_to_xdf_time(pupilData.fixations.(colName));
            
            for i = 1:height(pupilData.fixations)
                [timeDiff, sampleIdx] = min(abs(lslTimeInEEG - pupilData.fixations.timestamp_xdf(i)));
                
                if timeDiff < 0.1
                    currentEventCount = currentEventCount + 1;
                    amp.event(currentEventCount).latency = sampleIdx;
                    amp.event(currentEventCount).type = 'fixation';
                    amp.event(currentEventCount).source = 'pupil_fixations';
                    amp.event(currentEventCount).duration = 1;
                    
                    % Add fixation metadata
                    if ismember('durationms', pupilData.fixations.Properties.VariableNames)
                        duration_ms = pupilData.fixations.durationms(i);
                        if ~isnan(duration_ms)
                            amp.event(currentEventCount).duration = round(duration_ms * amp.srate / 1000);
                            amp.event(currentEventCount).fixation_duration_ms = duration_ms;
                        end
                    end
                    
                    % Add fixation position
                    if ismember('fixationXpx', pupilData.fixations.Properties.VariableNames)
                        amp.event(currentEventCount).fixation_x_px = pupilData.fixations.fixationXpx(i);
                    end
                    if ismember('fixationYpx', pupilData.fixations.Properties.VariableNames)
                        amp.event(currentEventCount).fixation_y_px = pupilData.fixations.fixationYpx(i);
                    end
                    
                    fixationsAdded = fixationsAdded + 1;
                end
            end
        end
    end
    
    totalEventsAdded = eventsAdded + blinksAdded + fixationsAdded;
    fprintf('   ✅ Added %d enhanced pupil events to EEG (Events:%d, Blinks:%d, Fixations:%d)\n', ...
           totalEventsAdded, eventsAdded, blinksAdded, fixationsAdded);
           
    %% Add enhanced gaze data to EEG channels from CSV
    fprintf('🔄 Adding enhanced gaze data from CSV to EEG channels...\n');
    
    if isfield(pupilData, 'gaze') && height(pupilData.gaze) > 0
        % Align gaze data timestamps
        gazeTimestampCol = find_timestamp_column(pupilData.gaze);
        if ~isempty(gazeTimestampCol)
            pupilData.gaze.timestamp_xdf = csv_to_xdf_time(pupilData.gaze.(gazeTimestampCol));
            
            % Find gaze X/Y and azimuth/elevation columns 
            gaze_x_col = '';
            gaze_y_col = '';
            azimuth_col = '';
            elevation_col = '';
            colNames = pupilData.gaze.Properties.VariableNames;
            
            % Look for gaze X column
            xCols = {'gazeXpx', 'gazeX_px_', 'gazeX', 'gaze_x'};
            for col = xCols
                if ismember(col{1}, colNames)
                    gaze_x_col = col{1};
                    break;
                end
            end
            
            % Look for gaze Y column
            yCols = {'gazeYpx', 'gazeY_px_', 'gazeY', 'gaze_y'};
            for col = yCols
                if ismember(col{1}, colNames)
                    gaze_y_col = col{1};
                    break;
                end
            end
            
            % Look for azimuth column
            azimuthCols = {'azimuth_deg_', 'azimuthdeg', 'azimuth', 'azimuth_deg'};
            for col = azimuthCols
                if ismember(col{1}, colNames)
                    azimuth_col = col{1};
                    break;
                end
            end
            
            % Look for elevation column
            elevationCols = {'elevation_deg_', 'elevationdeg', 'elevation', 'elevation_deg'};
            for col = elevationCols
                if ismember(col{1}, colNames)
                    elevation_col = col{1};
                    break;
                end
            end
            
            if ~isempty(gaze_x_col) && ~isempty(gaze_y_col)
                % Get EEG timing for interpolation
                lslTimeInEEG = amp.times / 1000; % Convert to seconds
                refIdx = find(amp.data(cntIdx,:) == streams{eegIdx}.time_series(1));
                lslTimeInEEG = lslTimeInEEG - lslTimeInEEG(refIdx) + streams{eegIdx}.time_stamps(1);
                
                % Filter gaze data to EEG time range
                eegTimeStart = min(lslTimeInEEG);
                eegTimeEnd = max(lslTimeInEEG);
                timeBuffer = 1.0; % 1 second buffer
                
                gazeTimeFilter = (pupilData.gaze.timestamp_xdf >= (eegTimeStart - timeBuffer)) & ...
                                (pupilData.gaze.timestamp_xdf <= (eegTimeEnd + timeBuffer));
                
                % Also remove invalid gaze points
                validGaze = gazeTimeFilter & ...
                           ~isnan(pupilData.gaze.(gaze_x_col)) & ...
                           ~isnan(pupilData.gaze.(gaze_y_col)) & ...
                           ~isnan(pupilData.gaze.timestamp_xdf);
                
                fprintf('   📊 Gaze samples: %d total, %d in time range, %d valid\n', ...
                       height(pupilData.gaze), sum(gazeTimeFilter), sum(validGaze));
                
                if sum(validGaze) > 10  % Need at least some valid data
                    % Extract filtered gaze data
                    gazeTimesValid = pupilData.gaze.timestamp_xdf(validGaze);
                    gazeXValid = pupilData.gaze.(gaze_x_col)(validGaze);
                    gazeYValid = pupilData.gaze.(gaze_y_col)(validGaze);
                    
                    % Interpolate gaze data to EEG counter times
                    gazeX = interp1(gazeTimesValid, gazeXValid, lslTimeInEEG, 'linear', NaN);
                    gazeY = interp1(gazeTimesValid, gazeYValid, lslTimeInEEG, 'linear', NaN);
                    
                    % Check interpolation quality
                    validInterp = ~isnan(gazeX) & ~isnan(gazeY);
                    interpCoverage = sum(validInterp) / length(validInterp) * 100;
                    
                    fprintf('   📊 Interpolation coverage: %.1f%% (%d/%d samples)\n', ...
                           interpCoverage, sum(validInterp), length(validInterp));
                    
                    % Determine how many channels to add
                    channelsToAdd = 2; % Start with X, Y
                    channelData = [gazeX; gazeY];
                    channelLabels = {'GazeX_CSV', 'GazeY_CSV'};
                    
                    % Add azimuth and elevation if available
                    if ~isempty(azimuth_col) && ismember(azimuth_col, colNames)
                        azimuthValid = pupilData.gaze.(azimuth_col)(validGaze);
                        azimuth = interp1(gazeTimesValid, azimuthValid, lslTimeInEEG, 'linear', NaN);
                        channelData = [channelData; azimuth];
                        channelLabels{end+1} = 'GazeAzimuth_CSV';
                        channelsToAdd = channelsToAdd + 1;
                        fprintf('   ✅ Found azimuth data in column: %s\n', azimuth_col);
                    end
                    
                    if ~isempty(elevation_col) && ismember(elevation_col, colNames)
                        elevationValid = pupilData.gaze.(elevation_col)(validGaze);
                        elevation = interp1(gazeTimesValid, elevationValid, lslTimeInEEG, 'linear', NaN);
                        channelData = [channelData; elevation];
                        channelLabels{end+1} = 'GazeElevation_CSV';
                        channelsToAdd = channelsToAdd + 1;
                        fprintf('   ✅ Found elevation data in column: %s\n', elevation_col);
                    end
                    
                    if interpCoverage > 1.0  % At least 1% coverage
                        % Add channels to EEG
                        newIdx = amp.nbchan + 1;
                        amp.data(newIdx:newIdx+channelsToAdd-1, :) = NaN(channelsToAdd, size(amp.data, 2));
                        amp.data(newIdx:newIdx+channelsToAdd-1, :) = channelData;
                        
                        % Set channel labels and types
                        for ch = 1:channelsToAdd
                            amp.chanlocs(newIdx + ch - 1).labels = channelLabels{ch};
                            amp.chanlocs(newIdx + ch - 1).type = 'Gaze';
                        end
                        
                        amp.nbchan = amp.nbchan + channelsToAdd;
                        
                        fprintf('   ✅ Added %d enhanced gaze channels to EEG (%.1f%% coverage)\n', channelsToAdd, interpCoverage);
                        
                        % Store gaze data for visualization
                        gazeDataForViz = struct();
                        gazeDataForViz.x = gazeXValid;
                        gazeDataForViz.y = gazeYValid;
                        gazeDataForViz.times = gazeTimesValid;
                        
                        % Store azimuth/elevation for potential future visualization
                        if exist('azimuthValid', 'var')
                            gazeDataForViz.azimuth = azimuthValid;
                        end
                        if exist('elevationValid', 'var')
                            gazeDataForViz.elevation = elevationValid;
                        end
                    else
                        fprintf('   ⚠️ Poor interpolation coverage (%.1f%%), skipping gaze channels\n', interpCoverage);
                    end
                else
                    fprintf('   ⚠️ Not enough valid gaze data in EEG time range (%d samples)\n', sum(validGaze));
                end
            else
                fprintf('   ⚠️ Could not find gaze X/Y columns. Available: %s\n', strjoin(colNames, ', '));
            end
        else
            fprintf('   ⚠️ Could not find timestamp column in gaze data\n');
        end
    else
        fprintf('   ℹ️ No gaze data available for enhanced processing\n');
    end
    
else
    fprintf('   ❌ Need at least 2 matched sync events for alignment\n');
end

%% After all data is merged and before returning, detect gait events
% Gait event detection using optimized algorithm
if ~isempty(xsensData)
    try
        fprintf('🚶 Detecting gait events from XSens foot data...\n');
        [heelStrikes, toeOffs, metrics] = gait_detection_final(amp, 'left', 'Verbose', true);
        
        % Convert gait events to EEGLAB event structure
        if ~isempty(heelStrikes) || ~isempty(toeOffs)
            newEvents = [];
            eventCount = 0;
            
            % Add heel strike events
            for i = 1:length(heelStrikes)
                eventCount = eventCount + 1;
                newEvents(eventCount).latency = round(heelStrikes(i) * amp.srate) + 1; % Convert to sample index
                newEvents(eventCount).duration = 0;
                newEvents(eventCount).type = 'heel_strike';
                newEvents(eventCount).source = 'xsens_gait_detection';
                newEvents(eventCount).pupil_type = '';
                newEvents(eventCount).pupil_timestamp_ns = NaN;
                newEvents(eventCount).pupil_timestamp_xdf = NaN;
                newEvents(eventCount).pupil_recording_id = '';
                newEvents(eventCount).blink_duration_ms = NaN;
                %newEvents(eventCount).urevent = eventCount;
            end
            
            % Add toe-off events
            for i = 1:length(toeOffs)
                eventCount = eventCount + 1;
                newEvents(eventCount).latency = round(toeOffs(i) * amp.srate) + 1; % Convert to sample index
                newEvents(eventCount).duration = 0;
                newEvents(eventCount).type = 'toe_off';
                newEvents(eventCount).source = 'xsens_gait_detection';
                newEvents(eventCount).pupil_type = '';
                newEvents(eventCount).pupil_timestamp_ns = NaN;
                newEvents(eventCount).pupil_timestamp_xdf = NaN;
                newEvents(eventCount).pupil_recording_id = '';
                newEvents(eventCount).blink_duration_ms = NaN;
                %newEvents(eventCount).urevent = eventCount;
            end
            
            % Merge with existing events
            if isfield(amp, 'event') && ~isempty(amp.event)
                originalEventCount = length(amp.event);
                for i = 1:length(newEvents)
                    %newEvents(i).urevent = originalEventCount + i;
                end
                amp.event = [amp.event, newEvents];
            else
                amp.event = newEvents;
            end
            
            fprintf('✅ Added %d gait events to EEG structure (%d heel strikes, %d toe-offs)\n', ...
                    eventCount, length(heelStrikes), length(toeOffs));
            
            % Store gait metrics in EEG structure for later analysis
            if ~isempty(fieldnames(metrics))
                amp.etc.gait_metrics = metrics;
                fprintf('📊 Gait metrics stored in EEG.etc.gait_metrics\n');
            end
        else
            fprintf('⚠️  No gait events detected\n');
        end
        
    catch ME
        fprintf('❌ Gait detection failed: %s\n', ME.message);
    end
end

%% Create comprehensive merging success visualizations
% re-order events
amp = eeg_checkset(amp,'eventconsistency');

% create merging metrics
fprintf('📊 Creating merging success visualizations...\n');

try
    % Create figure with multiple subplots
    fig = figure('Position', [100, 100, 1400, 1000], 'Name', sprintf('WaS4 Merge Summary - Subject %s', subjectFolder));
    
    % Subplot 1: Data timeline overview
    subplot(3, 3, 1);
    timelineData = [];
    timelineLabels = {};
    timelineColors = [];
    
    % EEG timeline
    if ~isempty(amp.times)
        eegTime = amp.times / 1000; % Convert to seconds
        timelineData(end+1,:) = [min(eegTime), max(eegTime)];
        timelineLabels{end+1} = 'EEG Data';
        timelineColors(end+1,:) = [0.2, 0.4, 0.8]; % Blue
    end
    
    % XSens timeline
    if ~isempty(xsensData) && any(xsensCoMIdx)
        xsensTime = streamsMerged{xsensCoMIdx}.time_stamps;
        timelineData(end+1,:) = [min(xsensTime), max(xsensTime)];
        timelineLabels{end+1} = 'XSens Data';
        timelineColors(end+1,:) = [0.8, 0.4, 0.2]; % Orange
    end
    
    % Pupil timeline
    if any(plGazeIdx)
        pupilTime = streamsMerged{plGazeIdx}.time_stamps;
        timelineData(end+1,:) = [min(pupilTime), max(pupilTime)];
        timelineLabels{end+1} = 'Pupil Data';
        timelineColors(end+1,:) = [0.2, 0.8, 0.4]; % Green
    end
    
    % ECG timeline
    if any(farosECGIdx)
        ecgTime = streamsMerged{farosECGIdx}.time_stamps;
        timelineData(end+1,:) = [min(ecgTime), max(ecgTime)];
        timelineLabels{end+1} = 'ECG Data';
        timelineColors(end+1,:) = [0.8, 0.2, 0.4]; % Red
    end
    
    % Plot timeline
    for i = 1:size(timelineData, 1)
        barh(i, timelineData(i,2) - timelineData(i,1), 'BaseValue', timelineData(i,1), ...
             'FaceColor', timelineColors(i,:), 'FaceAlpha', 0.7);
        hold on;
    end
    
    set(gca, 'YTick', 1:length(timelineLabels), 'YTickLabel', timelineLabels);
    xlabel('Time (seconds)');
    title('Data Stream Timeline');
    grid on;
    
    % Subplot 2: Channel count summary
    subplot(3, 3, 2);
    channelTypes = {};
    channelCounts = [];
    
    % Count channels by type
    eegChannels = sum(strcmp({amp.chanlocs.type}, 'EEG'));
    ecgChannels = sum(strcmp({amp.chanlocs.type}, 'ECG'));
    gazeChannels = sum(strcmp({amp.chanlocs.type}, 'Gaze'));
    xsensChannels = sum(strcmp({amp.chanlocs.type}, 'Xsens'));
    otherChannels = amp.nbchan - eegChannels - ecgChannels - gazeChannels - xsensChannels;
    
    channelTypes = {'EEG', 'ECG', 'Gaze', 'XSens', 'Other'};
    channelCounts = [eegChannels, ecgChannels, gazeChannels, xsensChannels, otherChannels];
    
    pie(channelCounts, channelTypes);
    title(sprintf('Channel Distribution (Total: %d)', amp.nbchan));
    
    % Subplot 3: Event summary
    subplot(3, 3, 3);
    if isfield(amp, 'event') && ~isempty(amp.event)
        eventSources = {};
        eventCounts = [];
        
        % Count events by source
        sources = {'pupil_csv_enhanced', 'gait_events', 'original'};
        for src = sources
            if isfield(amp.event, 'source')
                count = sum(strcmp({amp.event.source}, src{1}));
            else
                count = 0;
            end
            if count > 0
                eventSources{end+1} = src{1};
                eventCounts(end+1) = count;
            end
        end
        
        % Count events without source field (original events)
        noSourceCount = sum(~isfield(amp.event, 'source') | cellfun(@isempty, {amp.event.source}));
        if noSourceCount > 0
            eventSources{end+1} = 'original_events';
            eventCounts(end+1) = noSourceCount;
        end
        
        if ~isempty(eventCounts)
            bar(eventCounts);
            set(gca, 'XTickLabel', eventSources);
            xtickangle(45);
            title(sprintf('Event Count by Source (Total: %d)', length(amp.event)));
            ylabel('Number of Events');
        else
            text(0.5, 0.5, 'No events found', 'HorizontalAlignment', 'center');
            title('Event Summary');
        end
    else
        text(0.5, 0.5, 'No events found', 'HorizontalAlignment', 'center');
        title('Event Summary');
    end
    
    % Subplot 4: XSens data quality
    subplot(3, 3, 4);
    if ~isempty(xsensData) && any(xsensCoMIdx)
        comData = streamsMerged{xsensCoMIdx}.time_series;
        if size(comData, 1) >= 3
            plot(streamsMerged{xsensCoMIdx}.time_stamps, comData(1:3,:)');
            legend({'CoM X', 'CoM Y', 'CoM Z'}, 'Location', 'best');
            title('Center of Mass Position');
            xlabel('Time (s)');
            ylabel('Position (m)');
            grid on;
        end
    else
        text(0.5, 0.5, 'No XSens data available', 'HorizontalAlignment', 'center');
        title('XSens Data Quality');
    end
    
    % Subplot 5: ECG signal quality
    subplot(3, 3, 5);
    if any(farosECGIdx) && ~isempty(streamsMerged{farosECGIdx}.time_series)
        ecgSignal = streamsMerged{farosECGIdx}.time_series;
        ecgTime = streamsMerged{farosECGIdx}.time_stamps;
        
        % Show first 10 seconds of ECG
        maxTime = min(10, max(ecgTime) - min(ecgTime));
        timeIdx = ecgTime <= (min(ecgTime) + maxTime);
        
        plot(ecgTime(timeIdx), ecgSignal(timeIdx));
        title('ECG Signal Quality (First 10s)');
        xlabel('Time (s)');
        ylabel('Amplitude');
        grid on;
    else
        text(0.5, 0.5, 'No ECG data available', 'HorizontalAlignment', 'center');
        title('ECG Signal Quality');
    end
    
    % Subplot 6: Pupil gaze coverage
    subplot(3, 3, 6);
    if exist('gazeDataForViz', 'var') && isfield(gazeDataForViz, 'x') && isfield(gazeDataForViz, 'y')
        gazeX = gazeDataForViz.x;
        gazeY = gazeDataForViz.y;
        
        % Remove NaN values for plotting
        validIdx = ~isnan(gazeX) & ~isnan(gazeY);
        if sum(validIdx) > 0
            scatter(gazeX(validIdx), gazeY(validIdx), 1, 'filled', 'MarkerFaceAlpha', 0.3);
            title('Gaze Position Coverage (CSV)');
            xlabel('Gaze X (px)');
            ylabel('Gaze Y (px)');
            axis equal;
            grid on;
            
            % Add stats
            text(0.02, 0.98, sprintf('Valid samples: %d', sum(validIdx)), ...
                 'Units', 'normalized', 'VerticalAlignment', 'top', ...
                 'BackgroundColor', 'white', 'FontSize', 8);
        else
            text(0.5, 0.5, 'No valid gaze data', 'HorizontalAlignment', 'center');
            title('Gaze Position Coverage');
        end
    else
        text(0.5, 0.5, 'No gaze data available', 'HorizontalAlignment', 'center');
        title('Gaze Position Coverage');
    end
    
    % Subplot 7: Data synchronization quality
    subplot(3, 3, 7);
    syncQualityText = {};
    
    % Show sync information if available
    if exist('rmse', 'var')
        syncQualityText{end+1} = sprintf('Pupil Sync RMSE: %.6f s', rmse);
    end
    if exist('csvSyncEvents', 'var')
        syncQualityText{end+1} = sprintf('Pupil Sync Events: %d', length(csvSyncEvents));
    end
    
    if ~isempty(xsensData)
        syncQualityText{end+1} = 'XSens: Cross-correlation aligned';
    end
    
    syncQualityText{end+1} = sprintf('Total Channels: %d', amp.nbchan);
    syncQualityText{end+1} = sprintf('Total Events: %d', length(amp.event));
    syncQualityText{end+1} = sprintf('Recording Duration: %.1f min', (max(amp.times) - min(amp.times))/60000);
    
    text(0.1, 0.9, syncQualityText, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
         'FontSize', 10);
    title('Synchronization Summary');
    axis off;
    
    % Subplot 8: Data availability and gait detection overview
    subplot(3, 3, 8);
    
    % Check for missing data streams
    hasEEG = ~isempty(amp.data) && size(amp.data, 1) > 0;
    hasXSens = ~isempty(xsensData);
    hasPupil = any(plGazeIdx);
    hasECG = any(farosECGIdx);
    
    dataStatus = {'EEG', 'XSens', 'Pupil', 'ECG'};
    dataAvailable = [hasEEG, hasXSens, hasPupil, hasECG];
    
    colors = {'red', 'green'}; % red for missing, green for available
    statusText = {'✗ Missing', '✓ Available'}; % order matches boolean indexing
    
    % Data availability status (top half)
    for i = 1:length(dataStatus)
        colorIdx = dataAvailable(i) + 1; % 1 for missing (red), 2 for available (green)
        statusIdx = dataAvailable(i) + 1; % 1 for missing, 2 for available
        text(0.1, 0.95 - (i-1)*0.12, sprintf('%s: %s', dataStatus{i}, statusText{statusIdx}), ...
             'Units', 'normalized', 'Color', colors{colorIdx}, 'FontSize', 10, 'FontWeight', 'bold');
    end
    
    % Gait detection status (bottom half)
    text(0.1, 0.42, 'Gait Detection:', 'Units', 'normalized', 'FontSize', 11, 'FontWeight', 'bold');
    
    % Check for gait events in amp.event
    gaitEventsDetected = false;
    nHeelStrikes = 0;
    nToeOffs = 0;
    if isfield(amp, 'event') && ~isempty(amp.event)
        % Count gait events by type
        if isfield(amp.event, 'type')
            nHeelStrikes = sum(strcmp({amp.event.type}, 'heel_strike'));
            nToeOffs = sum(strcmp({amp.event.type}, 'toe_off'));
            gaitEventsDetected = (nHeelStrikes > 0) || (nToeOffs > 0);
        end
    end
    
    if hasXSens && gaitEventsDetected
        text(0.1, 0.32, sprintf('✓ %d heel strikes', nHeelStrikes), ...
             'Units', 'normalized', 'Color', 'green', 'FontSize', 9);
        text(0.1, 0.24, sprintf('✓ %d toe-offs', nToeOffs), ...
             'Units', 'normalized', 'Color', 'green', 'FontSize', 9);
             
        % Add gait metrics if available
        if isfield(amp, 'etc') && isfield(amp.etc, 'gait_metrics') && isfield(amp.etc.gait_metrics, 'cadence')
            text(0.1, 0.16, sprintf('Cadence: %.1f steps/min', amp.etc.gait_metrics.cadence), ...
                 'Units', 'normalized', 'Color', 'blue', 'FontSize', 9);
        end
    elseif hasXSens && ~gaitEventsDetected
        text(0.1, 0.32, '⚠ XSens available, no gait events', ...
             'Units', 'normalized', 'Color', [1 0.6 0], 'FontSize', 9); % Orange as RGB
    else
        text(0.1, 0.32, '✗ No XSens data for gait detection', ...
             'Units', 'normalized', 'Color', 'red', 'FontSize', 9);
    end
    
    title('Data Availability & Gait Detection');
    axis off;
    
    % Subplot 9: Processing summary
    subplot(3, 3, 9);
    summaryText = {
        sprintf('Subject: %s', subjectFolder),
        sprintf('Processing Date: %s', datestr(now)),
        sprintf('MATLAB Version: %s', version),
        sprintf('Function: %s', mfilename),
        '',
        'Data Integration Status:',
        sprintf('  • EEG channels merged: %d', eegChannels),
        sprintf('  • Physiological signals: %d', ecgChannels),
        sprintf('  • Gaze tracking: %d', gazeChannels),
        sprintf('  • Motion capture: %d', xsensChannels),
        sprintf('  • Total events: %d', length(amp.event))
    };
    
    text(0.05, 0.95, summaryText, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
         'FontSize', 9);
    title('Processing Summary');
    axis off;
    
    % Adjust layout
    sgtitle(sprintf('WaS4 Data Integration Summary - %s', subjectFolder), 'FontSize', 16, 'FontWeight', 'bold');
    
    % Save the figure
    figPath = fullfile(outPath, sprintf('%s_merge_summary.png', subjectFolder));
    saveas(fig, figPath);
    fprintf('   ✅ Saved visualization: %s\n', figPath);
    
    % Close figure to save memory
    close(fig);
    
catch ME
    warning('❌ Could not create visualizations: %s', ME.message);
end

fprintf('✅ Merging visualization complete!\n');

%% save data
eegMerged = amp;
eegMerged.nbchan = size(eegMerged.data,1);
eegMerged = eeg_checkset(eegMerged);
pop_saveset(eegMerged, 'filename', [subjectFolder '_merged.set'], 'filepath', fullfile(outPath));

end




function timealignedgaze = importGaze(filename)
%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 13);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["VarName1", "sectionId", "recordingId", "timestampns", "gazeXpx", "gazeYpx", "worn", "fixationId", "blinkId", "azimuthdeg", "elevationdeg", "timestamps", "lsl_times"];
opts.VariableTypes = ["double", "categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
% opts = setvaropts(opts, "blinkId", "WhitespaceRule", "preserve");
% opts = setvaropts(opts, ["sectionId", "recordingId", "blinkId"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["sectionId", "recordingId"], "EmptyFieldRule", "auto");

% Import the data
timealignedgaze = readtable(filename, opts);
end


function timestampCol = find_timestamp_column(dataTable)
    % Find the timestamp column in pupil labs CSV data
    timestampCol = '';
    
    % Possible timestamp column names in order of preference
    possibleCols = {'timestamp_ns_', 'timestamp_ns', 'timestamp', 'time', 'startTimestamp_ns_', 'startTimestampns'};
    
    for col = possibleCols
        if ismember(col{1}, dataTable.Properties.VariableNames)
            timestampCol = col{1};
            return;
        end
    end
    
    % If no standard column found, look for any column with 'timestamp' in the name
    colNames = dataTable.Properties.VariableNames;
    timestampCols = colNames(contains(colNames, 'timestamp', 'IgnoreCase', true));
    
    if ~isempty(timestampCols)
        timestampCol = timestampCols{1};
        return;
    end
    
    warning('Could not find timestamp column in data table');
end

function fixations = importFixations(filename)
%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 10);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["sectionId", "recordingId", "fixationId", "startTimestampns", "endTimestampns", "durationms", "fixationXpx", "fixationYpx", "azimuthdeg", "elevationdeg"];
opts.VariableTypes = ["categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["sectionId", "recordingId"], "EmptyFieldRule", "auto");

% Import the data
fixations = readtable(filename, opts);
end

%% Time alignment function for Pupil Labs data
function success = create_time_aligned_gaze(pupil_folder_path)
% CREATE_TIME_ALIGNED_GAZE Create time_aligned_gaze.csv from gaze.csv and time_alignment_parameters.json
%
% This function performs the same time alignment as the Python script:
% 1. Load gaze.csv and time_alignment_parameters.json
% 2. Convert timestamps from nanoseconds to seconds
% 3. Apply linear mapping using intercept and slope from alignment JSON
% 4. Save the result as time_aligned_gaze.csv

success = false;

try
    % File paths
    gaze_file = fullfile(pupil_folder_path, 'gaze.csv');
    params_file = fullfile(pupil_folder_path, 'time_alignment_parameters.json');
    output_file = fullfile(pupil_folder_path, 'time_aligned_gaze.csv');
    
    % Check if required files exist
    if ~exist(gaze_file, 'file')
        fprintf('❌ gaze.csv not found in %s\n', pupil_folder_path);
        return;
    end
    
    if ~exist(params_file, 'file')
        fprintf('❌ time_alignment_parameters.json not found in %s\n', pupil_folder_path);
        return;
    end
    
    fprintf('🔗 Creating time-aligned gaze data from %s\n', pupil_folder_path);
    
    % Load gaze data
    fprintf('   Loading gaze.csv...\n');
    gaze_data = readtable(gaze_file);
    
    % Convert nanoseconds to seconds (same as Python: * 1e-9)
    if ismember('timestamp__ns_', gaze_data.Properties.VariableNames)
        timestamp_col = 'timestamp__ns_';
    elseif ismember('timestamp [ns]', gaze_data.Properties.VariableNames)
        timestamp_col = 'timestamp [ns]';
    else
        fprintf('❌ Could not find timestamp column in gaze.csv\n');
        return;
    end
    
    gaze_data.timestamp_s = gaze_data.(timestamp_col) * 1e-9;
    
    % Load alignment parameters
    fprintf('   Loading time_alignment_parameters.json...\n');
    json_text = fileread(params_file);
    params = jsondecode(json_text);
    
    % Apply linear mapping: lsl_time = intercept + timestamp_s * slope
    % This matches the Python code: perform_linear_mapping(cloud_gaze_data["timestamp [s]"], parameter_dict["cloud_to_lsl"])
    intercept = params.cloud_to_lsl.intercept;
    slope = params.cloud_to_lsl.slope;
    
    gaze_data.lsl_time_s = intercept + gaze_data.timestamp_s * slope;
    
    fprintf('   Applied linear mapping: intercept=%.6f, slope=%.6f\n', intercept, slope);
    
    % Save the time-aligned data
    fprintf('   Saving time_aligned_gaze.csv...\n');
    writetable(gaze_data, output_file);
    
    fprintf('✅ Successfully created %s\n', output_file);
    success = true;
    
catch ME
    fprintf('❌ Time alignment failed: %s\n', ME.message);
    success = false;
end

end
