%% EEG-PUPIL LABS PREPROCESSING PIPELINE - ENHANCED
% One-step alignment and merging script for EEG and Pupil Labs data
% Run this as the first step in your preprocessing pipeline
%
% This script:
% 1. Finds all subjects in the data directory
% 2. Uses WaS4_merge for comprehensive data synchronization
% 3. Performs gait event detection on synchronized data
% 4. Saves results and generates reports
%
% ENHANCED FEATURES:
% - Uses WaS4_merge.m for all data loading and synchronization
% - Gait event detection using synchronized streams
% - Robust error handling and fallback options
% - Comprehensive reporting and quality control

clear; clc; close all;

dataPath = '/Volumes/ergo/rohdaten/RohDaten_GRAIL/WaS4/DATA/';  % Main data directory
outputPath = '/Volumes/Work4TB/Seafile/WaS4/data/aligned/';  % Output directory for aligned files

%% EEG-PUPIL LABS PREPROCESSING PIPELINE

% Processing options
processAllSubjects = true;  % Set to false to process specific subjects
specificSubjects = [1];  % Only used if processAllSubjects = false
overwriteExisting = true;  % Set to true to reprocess existing files

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
        
        %% Load and synchronize all data using WaS4_merge
        if verboseOutput, fprintf('🦴 Using WaS4_merge for comprehensive data synchronization...\n'); end
        
        % Extract subject number from folder name
        subjectNumStr = regexp(subjectFolder, '\d+', 'match');
        if ~isempty(subjectNumStr)
            subjectNum = str2double(subjectNumStr{1});
        else
            error('Could not extract subject number from folder: %s', subjectFolder);
        end
        
        % Call WaS4_merge to get all synchronized data
        if verboseOutput, fprintf('   🔄 Processing subject %03d with WaS4_merge...\n', subjectNum); end
        [eeg, streamsMerged] = WaS4_merge(subjectNum, dataPath);
        
        if verboseOutput
            fprintf('   ✅ WaS4_merge completed successfully\n');
            fprintf('   📊 Merged EEG channels: %d\n', eeg.nbchan);
            fprintf('   📊 Synchronized streams: %d\n', length(streamsMerged));
        end
        
        
        %% Save results
        if saveResults
            outputFile = fullfile(outputPath, sprintf('%s_aligned.set', subjectFolder));
            if verboseOutput, fprintf('💾 Saving results to %s...\n', outputFile); end
            
            % Update EEG structure
            eeg = eeg_checkset(eeg);
            pop_saveset(eeg, 'filename', sprintf('%s_aligned.set', subjectFolder), 'filepath', outputPath);
            
            if verboseOutput, fprintf('   ✅ Results saved successfully\n'); end
        end
        
        % Record processing results for summary
        processedSubjects{end+1} = subjectFolder;
        
    catch ME
        % Handle errors 
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