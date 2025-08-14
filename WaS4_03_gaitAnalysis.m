%% GAIT ANALYSIS SCRIPT - WaS4 Study
% Detects gait events and performs gait analysis on merged EEG data
% Run this script AFTER WaS4_02_edfComMerge.m
%
% This script:
% 1. Loads merged EEG files with XSens foot position data
% 2. Detects heel-on, toe-on, heel-off, toe-off events
% 3. Calculates gait parameters (stride duration, cadence, etc.)
% 4. Adds gait events to EEG structure for further analysis
% 5. Saves enhanced EEG files with gait events
% 6. Creates comprehensive gait analysis plots
%
% DEPENDENCIES: WaS4_02_edfComMerge.m, detect_gait_events.m

clear; clc; close all;

%% ==================== CONFIGURATION ====================

% Data paths
mergedPath = '/Volumes/Work4TB/Seafile/WaS4/data/merged/';    % Merged data directory (from WaS4_02)
gaitPath = '/Volumes/Work4TB/Seafile/WaS4/data/gait/';        % Gait analysis output directory

% Processing options
processAllSubjects = true;  % Set to false to process specific subjects
specificSubjects = [46];    % Only used if processAllSubjects = false
overwriteExisting = false;  % Set to true to reprocess existing files

% Gait detection parameters
gaitParams = struct();
gaitParams.minStepDuration = 0.4;      % Minimum step duration in seconds
gaitParams.maxStepDuration = 2.0;      % Maximum step duration in seconds
gaitParams.heightThreshold = 0.02;     % Foot height threshold in meters (2cm)
gaitParams.velocityThreshold = 0.1;    % Velocity threshold for event detection
gaitParams.filterCutoff = 10;          % Low-pass filter cutoff in Hz
gaitParams.plotResults = true;         % Create diagnostic plots
gaitParams.verboseOutput = true;       % Detailed console output

% Output options
saveResults = true;
savePlots = true;
verboseOutput = true;

%% ==================== INITIALIZATION ====================

if verboseOutput
    fprintf('==========================================================\n');
    fprintf('GAIT ANALYSIS SCRIPT - WaS4 Study\n');
    fprintf('==========================================================\n\n');
    fprintf('📁 Merged data path: %s\n', mergedPath);
    fprintf('📁 Gait analysis output: %s\n', gaitPath);
    fprintf('🦶 Gait detection parameters:\n');
    fprintf('   - Step duration range: %.1f - %.1f seconds\n', ...
           gaitParams.minStepDuration, gaitParams.maxStepDuration);
    fprintf('   - Height threshold: %.2f meters\n', gaitParams.heightThreshold);
    fprintf('   - Filter cutoff: %.1f Hz\n', gaitParams.filterCutoff);
    fprintf('\n');
end

% Check directories
if ~exist(mergedPath, 'dir')
    error('Merged data directory not found: %s (run WaS4_02_edfComMerge.m first)', mergedPath);
end
if ~exist(gaitPath, 'dir')
    mkdir(gaitPath);
    if verboseOutput, fprintf('📁 Created gait analysis directory: %s\n\n', gaitPath); end
end

% Find merged EEG files
mergedFiles = dir(fullfile(mergedPath, '*_merged.set'));
if isempty(mergedFiles)
    error('No merged EEG files found in %s (run WaS4_02_edfComMerge.m first)', mergedPath);
end

% Extract subject list from merged files
subjectList = {};
for i = 1:length(mergedFiles)
    filename = mergedFiles(i).name;
    subjectMatch = regexp(filename, '(WaS_\d+)_merged\.set', 'tokens', 'once');
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
            warning('Subject %03d not found in merged data', subjectNum);
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
gaitAnalysisStats = [];

%% ==================== MAIN PROCESSING LOOP ====================

for s = 1:length(subjectList)
    subjectFolder = subjectList{s};
    
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
        fprintf('🚶 ANALYZING GAIT: %s (Subject %d) [%d/%d]\n', subjectFolder, subjectNum, s, length(subjectList));
        fprintf('==========================================================\n\n');
    end
    
    try
        % Check if already processed
        outputFile = fullfile(gaitPath, sprintf('%s_gait.set', subjectFolder));
        if exist(outputFile, 'file') && ~overwriteExisting
            if verboseOutput
                fprintf('⏭️  Already processed (use overwriteExisting=true to reprocess)\n\n');
            end
            skippedSubjects{end+1} = subjectFolder;
            continue;
        end
        
        %% Step 1: Load merged EEG data
        if verboseOutput, fprintf('🧠 Loading merged EEG data...\n'); end
        
        mergedFile = fullfile(mergedPath, sprintf('%s_merged.set', subjectFolder));
        if ~exist(mergedFile, 'file')
            error('Merged EEG file not found: %s', mergedFile);
        end
        
        eeg = pop_loadset('filename', sprintf('%s_merged.set', subjectFolder), 'filepath', mergedPath);
        
        if verboseOutput
            fprintf('   ✅ Loaded merged EEG: %d channels, %d samples, %.1f Hz\n', ...
                   eeg.nbchan, eeg.pnts, eeg.srate);
            
            % Show available channel types
            channelTypes = {eeg.chanlocs.type};
            uniqueTypes = unique(channelTypes);
            uniqueTypes = uniqueTypes(~cellfun(@isempty, uniqueTypes));
            fprintf('   📊 Available channel types: %s\n', strjoin(uniqueTypes, ', '));
        end
        
        %% Step 2: Check for required XSens foot position data
        channelLabels = {eeg.chanlocs.labels};
        footChannels = find(~cellfun(@isempty, regexpi(channelLabels, '.*foot.*')));
        
        if isempty(footChannels)
            error('No foot position channels found. Ensure XSens data is included in WaS4_02_edfComMerge.m');
        end
        
        if verboseOutput
            fprintf('   🦶 Found foot-related channels:\n');
            for f = 1:length(footChannels)
                chIdx = footChannels(f);
                fprintf('      - Ch %d: %s (%s)\n', chIdx, channelLabels{chIdx}, eeg.chanlocs(chIdx).type);
            end
            fprintf('\n');
        end
        
        %% Step 3: Detect gait events
        if verboseOutput, fprintf('🔄 Detecting gait events...\n'); end
        
        startTime = tic;
        [gaitEvents, gaitStats] = detect_gait_events(eeg, ...
            'minStepDuration', gaitParams.minStepDuration, ...
            'maxStepDuration', gaitParams.maxStepDuration, ...
            'heightThreshold', gaitParams.heightThreshold, ...
            'velocityThreshold', gaitParams.velocityThreshold, ...
            'filterCutoff', gaitParams.filterCutoff, ...
            'plotResults', false, ... % We'll create our own plots
            'verboseOutput', gaitParams.verboseOutput);
        
        processingTime = toc(startTime);
        
        if verboseOutput
            fprintf('   ✅ Gait event detection completed in %.1f seconds\n\n', processingTime);
        end
        
        %% Step 4: Add gait events to EEG structure
        if verboseOutput, fprintf('🔄 Adding gait events to EEG structure...\n'); end
        
        % Count existing events
        initialEventCount = length(eeg.event);
        
        % Add gait events for each foot
        feet = {'left', 'right'};
        eventTypes = {'heel_on', 'toe_on', 'heel_off', 'toe_off'};
        eventsAdded = 0;
        
        for f = 1:length(feet)
            foot = feet{f};
            if isfield(gaitEvents, foot)
                for e = 1:length(eventTypes)
                    eventType = eventTypes{e};
                    if isfield(gaitEvents.(foot), eventType) && ~isempty(gaitEvents.(foot).(eventType))
                        eventSamples = gaitEvents.(foot).(eventType);
                        
                        for i = 1:length(eventSamples)
                            eventIdx = length(eeg.event) + 1;
                            eeg.event(eventIdx).latency = eventSamples(i);
                            eeg.event(eventIdx).duration = 1;
                            eeg.event(eventIdx).type = sprintf('%s_%s', foot, eventType);
                            eeg.event(eventIdx).source = 'gait_analysis';
                            eeg.event(eventIdx).foot = foot;
                            eeg.event(eventIdx).gait_phase = eventType;
                            eventsAdded = eventsAdded + 1;
                        end
                    end
                end
            end
        end
        
        % Re-sort events by latency
        if eventsAdded > 0
            eeg = eeg_checkset(eeg, 'eventconsistency');
        end
        
        if verboseOutput
            fprintf('   ✅ Added %d gait events to EEG structure\n', eventsAdded);
            fprintf('   📊 Total events: %d → %d\n', initialEventCount, length(eeg.event));
        end
        
        %% Step 5: Store analysis statistics
        analysisStats = struct();
        analysisStats.subject = subjectFolder;
        analysisStats.processing_time_seconds = processingTime;
        analysisStats.events_added = eventsAdded;
        analysisStats.gait_stats = gaitStats;
        analysisStats.detection_parameters = gaitParams;
        analysisStats.analysis_date = datetime('now');
        
        gaitAnalysisStats = [gaitAnalysisStats; analysisStats];
        
        %% Step 6: Create comprehensive gait plots
        if gaitParams.plotResults
            if verboseOutput, fprintf('\n📈 Creating gait analysis plots...\n'); end
            
            try
                % Set up figure for this subject
                figTitle = sprintf('Gait Analysis - %s', subjectFolder);
                create_comprehensive_gait_plots(eeg, gaitEvents, gaitStats, subjectFolder, ...
                                               gaitPath, savePlots, verboseOutput);
                
                if verboseOutput
                    fprintf('   ✅ Gait plots created\n');
                end
            catch ME
                if verboseOutput
                    fprintf('   ⚠️ Error creating plots: %s\n', ME.message);
                end
            end
        end
        
        %% Step 7: Save results
        if saveResults
            if verboseOutput, fprintf('\n💾 Saving gait analysis results...\n'); end
            
            % Save enhanced EEG file with gait events
            pop_saveset(eeg, 'filename', sprintf('%s_gait.set', subjectFolder), 'filepath', gaitPath);
            
            % Save detailed gait analysis data
            gaitDataFile = fullfile(gaitPath, sprintf('%s_gait_analysis.mat', subjectFolder));
            gaitAnalysisData = struct();
            gaitAnalysisData.gaitEvents = gaitEvents;
            gaitAnalysisData.gaitStats = gaitStats;
            gaitAnalysisData.analysisStats = analysisStats;
            gaitAnalysisData.detectionParameters = gaitParams;
            save(gaitDataFile, 'gaitAnalysisData', '-v7.3');
            
            if verboseOutput
                fprintf('   ✅ Enhanced EEG saved: %s\n', outputFile);
                fprintf('   ✅ Gait analysis data saved: %s\n', gaitDataFile);
            end
        end
        
        % Mark as successfully processed
        processedSubjects{end+1} = subjectFolder;
        
        if verboseOutput
            fprintf('\n✅ %s GAIT ANALYSIS COMPLETED!\n', subjectFolder);
            fprintf('   🚶 Total steps detected: %d\n', gaitStats.total_steps);
            if gaitStats.total_steps > 0
                fprintf('   📊 Step frequency: %.1f steps/min\n', gaitStats.step_frequency);
                if isfield(gaitStats, 'stride_duration')
                    fprintf('   ⏱️  Mean stride duration: %.2f ± %.2f seconds\n', ...
                           gaitStats.stride_duration.mean, gaitStats.stride_duration.std);
                end
            end
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
    fprintf('GAIT ANALYSIS SUMMARY\n');
    fprintf('==========================================================\n\n');
    
    fprintf('✅ SUCCESSFULLY PROCESSED (%d subjects):\n', length(processedSubjects));
    for i = 1:length(processedSubjects)
        if i <= length(gaitAnalysisStats)
            stats = gaitAnalysisStats(i);
            fprintf('   %s (%d steps, %.1f/min, %.1fs processing)\n', processedSubjects{i}, ...
                   stats.gait_stats.total_steps, stats.gait_stats.step_frequency, ...
                   stats.processing_time_seconds);
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
    if ~isempty(gaitAnalysisStats)
        allSteps = [gaitAnalysisStats.gait_stats];
        totalSteps = sum([allSteps.total_steps]);
        validSubjects = sum([allSteps.total_steps] > 0);
        
        fprintf('\n📊 OVERALL GAIT STATISTICS:\n');
        fprintf('   Total steps detected: %d\n', totalSteps);
        fprintf('   Subjects with gait data: %d/%d (%.1f%%)\n', ...
               validSubjects, length(processedSubjects), ...
               validSubjects/length(processedSubjects)*100);
        
        if validSubjects > 0
            validStats = gaitAnalysisStats([gaitAnalysisStats.gait_stats.total_steps] > 0);
            stepFreqs = [validStats.gait_stats];
            stepFreqs = [stepFreqs.step_frequency];
            
            fprintf('   Average step frequency: %.1f ± %.1f steps/min\n', ...
                   mean(stepFreqs), std(stepFreqs));
            
            % Stride duration statistics (if available)
            hasStrideDuration = arrayfun(@(x) isfield(x.gait_stats, 'stride_duration'), validStats);
            if any(hasStrideDuration)
                strideStats = validStats(hasStrideDuration);
                strideMeans = arrayfun(@(x) x.gait_stats.stride_duration.mean, strideStats);
                fprintf('   Average stride duration: %.2f ± %.2f seconds\n', ...
                       mean(strideMeans), std(strideMeans));
            end
        end
        
        totalProcessingTime = sum([gaitAnalysisStats.processing_time_seconds]);
        fprintf('\n⏱️  PROCESSING PERFORMANCE:\n');
        fprintf('   Total processing time: %.1f seconds (%.1f minutes)\n', ...
               totalProcessingTime, totalProcessingTime/60);
        fprintf('   Average per subject: %.1f seconds\n', ...
               totalProcessingTime / length(processedSubjects));
    end
    
    fprintf('\n💾 Gait analysis output: %s\n', gaitPath);
    fprintf('==========================================================\n');
end

%% Helper function for comprehensive plotting
function create_comprehensive_gait_plots(eeg, gaitEvents, gaitStats, subjectName, outputPath, savePlots, verbose)
    try
        % Create main gait analysis figure
        fig = figure('Position', [50, 50, 1800, 1200], 'Name', sprintf('%s - Comprehensive Gait Analysis', subjectName));
        
        timeVector = eeg.times / 1000; % Convert to seconds
        feet = {'left', 'right'};
        colors = {[0.2, 0.4, 0.8], [0.8, 0.2, 0.4]}; % Blue for left, red for right
        
        % Plot 1: Foot trajectories with all gait events
        subplot(3, 3, 1);
        hold on;
        
        for f = 1:2
            foot = feet{f};
            if isfield(gaitEvents, foot) && isfield(gaitEvents.(foot), 'foot_height')
                % Plot foot height
                plot(timeVector, gaitEvents.(foot).foot_height, 'Color', colors{f}, ...
                     'LineWidth', 2, 'DisplayName', sprintf('%s foot', foot));
                
                % Mark all gait events
                eventTypes = {'heel_on', 'toe_on', 'heel_off', 'toe_off'};
                markers = {'^', 'o', 'v', 's'};
                eventColors = {'g', 'b', 'r', 'm'};
                
                for e = 1:length(eventTypes)
                    eventType = eventTypes{e};
                    if isfield(gaitEvents.(foot), eventType) && ~isempty(gaitEvents.(foot).(eventType))
                        eventSamples = gaitEvents.(foot).(eventType);
                        eventTimes = timeVector(eventSamples);
                        eventHeights = gaitEvents.(foot).foot_height(eventSamples);
                        
                        scatter(eventTimes, eventHeights, 80, eventColors{e}, markers{e}, ...
                               'filled', 'DisplayName', sprintf('%s %s', foot, strrep(eventType, '_', '-')));
                    end
                end
            end
        end
        
        title('Foot Height with Gait Events');
        xlabel('Time (s)');
        ylabel('Height (m)');
        legend('Location', 'eastoutside');
        grid on;
        
        % Plot 2: Step frequency over time
        subplot(3, 3, 2);
        if gaitStats.total_steps > 5
            % Calculate local step frequency using sliding window
            windowSize = 30; % seconds
            stepTimes = [];
            
            for f = 1:2
                foot = feet{f};
                if isfield(gaitEvents, foot) && isfield(gaitEvents.(foot), 'heel_on')
                    heelOnTimes = timeVector(gaitEvents.(foot).heel_on);
                    stepTimes = [stepTimes, heelOnTimes];
                end
            end
            
            if length(stepTimes) > 5
                stepTimes = sort(stepTimes);
                timePoints = min(stepTimes):5:max(stepTimes); % Every 5 seconds
                localFreq = zeros(size(timePoints));
                
                for t = 1:length(timePoints)
                    windowStart = timePoints(t) - windowSize/2;
                    windowEnd = timePoints(t) + windowSize/2;
                    stepsInWindow = sum(stepTimes >= windowStart & stepTimes <= windowEnd);
                    localFreq(t) = stepsInWindow / windowSize * 60; % steps/min
                end
                
                plot(timePoints, localFreq, 'k-', 'LineWidth', 2);
                title(sprintf('Step Frequency Over Time\n(Mean: %.1f steps/min)', gaitStats.step_frequency));
                xlabel('Time (s)');
                ylabel('Steps/min');
                grid on;
            else
                text(0.5, 0.5, 'Insufficient data for frequency analysis', ...
                     'Units', 'normalized', 'HorizontalAlignment', 'center');
                title('Step Frequency Over Time');
            end
        else
            text(0.5, 0.5, 'Insufficient steps detected', ...
                 'Units', 'normalized', 'HorizontalAlignment', 'center');
            title('Step Frequency Over Time');
        end
        
        % Plot 3: Stride duration analysis
        subplot(3, 3, 3);
        allCycles = [];
        for f = 1:2
            foot = feet{f};
            if isfield(gaitEvents, foot) && isfield(gaitEvents.(foot), 'step_cycles')
                allCycles = [allCycles; gaitEvents.(foot).step_cycles];
            end
        end
        
        if ~isempty(allCycles) && length(allCycles) > 2
            strideDurations = [allCycles.stride_duration];
            histogram(strideDurations, max(5, round(length(strideDurations)/3)), ...
                     'FaceColor', [0.6, 0.6, 0.8], 'EdgeColor', 'k', 'FaceAlpha', 0.7);
            title(sprintf('Stride Duration\n%.2f±%.2f s (CV: %.1f%%)', ...
                         mean(strideDurations), std(strideDurations), ...
                         std(strideDurations)/mean(strideDurations)*100));
            xlabel('Duration (s)');
            ylabel('Count');
            grid on;
        else
            text(0.5, 0.5, 'No stride data available', ...
                 'Units', 'normalized', 'HorizontalAlignment', 'center');
            title('Stride Duration Distribution');
        end
        
        % Plot 4-5: Individual foot gait phases
        for f = 1:2
            subplot(3, 3, 3+f);
            foot = feet{f};
            
            if isfield(gaitEvents, foot) && isfield(gaitEvents.(foot), 'step_cycles')
                cycles = gaitEvents.(foot).step_cycles;
                
                if ~isempty(cycles) && length(cycles) > 1
                    % Create phase diagram
                    for i = 1:min(length(cycles), 20) % Limit to first 20 cycles for clarity
                        c = cycles(i);
                        
                        % Stance phase (heel-on to toe-off)
                        stanceWidth = c.toe_off_time - c.heel_on_time;
                        
                        % Only plot rectangle if width is positive
                        if stanceWidth > 0 && ~isnan(stanceWidth) && ~isnan(c.heel_on_time)
                            rectangle('Position', [c.heel_on_time, i-0.4, stanceWidth, 0.8], ...
                                     'FaceColor', colors{f}, 'EdgeColor', 'k', 'FaceAlpha', 0.6);
                        end
                        
                        % Mark heel-off and toe-on if available
                        if isfield(c, 'heel_off_time') && ~isnan(c.heel_off_time)
                            line([c.heel_off_time, c.heel_off_time], [i-0.4, i+0.4], ...
                                 'Color', 'w', 'LineWidth', 3);
                        end
                        
                        if isfield(c, 'toe_on_time') && ~isnan(c.toe_on_time)
                            line([c.toe_on_time, c.toe_on_time], [i-0.4, i+0.4], ...
                                 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
                        end
                    end
                    
                    ylim([0.5, min(length(cycles), 20)+0.5]);
                    title(sprintf('%s Foot Gait Cycles (n=%d)', ...
                                 upper(foot(1))+lower(foot(2:end)), length(cycles)));
                    xlabel('Time (s)');
                    ylabel('Cycle #');
                    grid on;
                else
                    text(0.5, 0.5, sprintf('No %s foot cycles', foot), ...
                         'Units', 'normalized', 'HorizontalAlignment', 'center');
                    title(sprintf('%s Foot Gait Cycles', upper(foot(1))+lower(foot(2:end))));
                end
            else
                text(0.5, 0.5, sprintf('No %s foot data', foot), ...
                     'Units', 'normalized', 'HorizontalAlignment', 'center');
                title(sprintf('%s Foot Gait Cycles', upper(foot(1))+lower(foot(2:end))));
            end
        end
        
        % Plot 6: Stance vs Swing phase analysis
        subplot(3, 3, 6);
        if ~isempty(allCycles)
            stanceDurations = [allCycles.stance_duration];
            swingDurations = [allCycles.swing_duration];
            
            scatter(stanceDurations, swingDurations, 50, [0.5, 0.5, 0.5], 'filled', 'MarkerFaceAlpha', 0.6);
            hold on;
            
            % Add regression line
            if length(stanceDurations) > 3
                p = polyfit(stanceDurations, swingDurations, 1);
                xfit = linspace(min(stanceDurations), max(stanceDurations), 100);
                yfit = polyval(p, xfit);
                plot(xfit, yfit, 'r-', 'LineWidth', 2);
            end
            
            xlabel('Stance Duration (s)');
            ylabel('Swing Duration (s)');
            title(sprintf('Stance vs Swing\nStance: %.1f%%, Swing: %.1f%%', ...
                         gaitStats.stance_percent, gaitStats.swing_percent));
            grid on;
            axis equal;
        else
            text(0.5, 0.5, 'No phase data available', ...
                 'Units', 'normalized', 'HorizontalAlignment', 'center');
            title('Stance vs Swing Phase');
        end
        
        % Plot 7: Event timeline
        subplot(3, 3, 7);
        yPos = 0;
        legendEntries = {};
        
        for f = 1:2
            foot = feet{f};
            if isfield(gaitEvents, foot)
                eventTypes = {'heel_on', 'toe_on', 'heel_off', 'toe_off'};
                eventColors = {'g', 'b', 'r', 'm'};
                eventMarkers = {'^', 'o', 'v', 's'};
                
                for e = 1:length(eventTypes)
                    eventType = eventTypes{e};
                    if isfield(gaitEvents.(foot), eventType) && ~isempty(gaitEvents.(foot).(eventType))
                        eventTimes = timeVector(gaitEvents.(foot).(eventType));
                        scatter(eventTimes, repmat(yPos, size(eventTimes)), 60, ...
                               eventColors{e}, eventMarkers{e}, 'filled');
                        hold on;
                        
                        if f == 1 % Only add to legend once
                            legendEntries{end+1} = strrep(eventType, '_', '-');
                        end
                    end
                    yPos = yPos + 1;
                end
                yPos = yPos + 1; % Gap between feet
            end
        end
        
        title('Gait Event Timeline');
        xlabel('Time (s)');
        ylabel('Event Type');
        if ~isempty(legendEntries)
            legend(legendEntries, 'Location', 'eastoutside');
        end
        grid on;
        
        % Plot 8: Summary statistics
        subplot(3, 3, 8);
        axis off;
        
        summaryText = {
            sprintf('GAIT SUMMARY: %s', subjectName),
            '',
            sprintf('Total Steps: %d', gaitStats.total_steps),
            sprintf('Step Frequency: %.1f steps/min', gaitStats.step_frequency),
            ''
        };
        
        if isfield(gaitStats, 'stride_duration') && gaitStats.total_steps > 0
            summaryText{end+1} = 'TEMPORAL PARAMETERS:';
            summaryText{end+1} = sprintf('Stride: %.2f±%.2f s', ...
                                       gaitStats.stride_duration.mean, gaitStats.stride_duration.std);
            summaryText{end+1} = sprintf('Variability: %.1f%%', ...
                                       gaitStats.stride_duration.cv * 100);
            
            if isfield(gaitStats, 'stance_percent')
                summaryText{end+1} = '';
                summaryText{end+1} = 'PHASE ANALYSIS:';
                summaryText{end+1} = sprintf('Stance: %.1f%%', gaitStats.stance_percent);
                summaryText{end+1} = sprintf('Swing: %.1f%%', gaitStats.swing_percent);
            end
        end
        
        % Add foot-specific information
        for f = 1:2
            foot = feet{f};
            if isfield(gaitEvents, foot)
                summaryText{end+1} = '';
                summaryText{end+1} = sprintf('%s FOOT:', upper(foot));
                
                eventTypes = {'heel_on', 'toe_on', 'heel_off', 'toe_off'};
                for e = 1:length(eventTypes)
                    eventType = eventTypes{e};
                    if isfield(gaitEvents.(foot), eventType)
                        count = length(gaitEvents.(foot).(eventType));
                        summaryText{end+1} = sprintf('  %s: %d', strrep(eventType, '_', '-'), count);
                    end
                end
            end
        end
        
        text(0.05, 0.95, summaryText, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
             'FontSize', 9, 'FontWeight', 'normal');
        
        % Plot 9: Data quality indicators
        subplot(3, 3, 9);
        qualityMetrics = [];
        qualityLabels = {};
        
        % Calculate data coverage for each foot
        for f = 1:2
            foot = feet{f};
            if isfield(gaitEvents, foot) && isfield(gaitEvents.(foot), 'foot_height')
                footData = gaitEvents.(foot).foot_height;
                coverage = sum(~isnan(footData)) / length(footData) * 100;
                qualityMetrics(end+1) = coverage;
                qualityLabels{end+1} = sprintf('%s foot coverage', foot);
            end
        end
        
        % Add step detection success rate
        if gaitStats.total_steps > 0
            recordingDuration = timeVector(end) - timeVector(1);
            expectedSteps = recordingDuration / 60 * 120; % Assume ~120 steps/min normal walking
            detectionRate = min(100, gaitStats.total_steps / expectedSteps * 100);
            qualityMetrics(end+1) = detectionRate;
            qualityLabels{end+1} = 'Step detection rate';
        end
        
        if ~isempty(qualityMetrics)
            bar(qualityMetrics);
            set(gca, 'XTickLabel', qualityLabels, 'XTickLabelRotation', 45);
            title('Data Quality Metrics');
            ylabel('Percentage (%)');
            ylim([0, 100]);
            grid on;
        else
            text(0.5, 0.5, 'No quality metrics available', ...
                 'Units', 'normalized', 'HorizontalAlignment', 'center');
            title('Data Quality Metrics');
        end
        
        sgtitle(sprintf('Comprehensive Gait Analysis - %s', subjectName), 'FontSize', 14, 'FontWeight', 'bold');
        
        if savePlots
            plotPath = fullfile(outputPath, sprintf('%s_comprehensive_gait_analysis.png', subjectName));
            print(fig, plotPath, '-dpng', '-r300');
            if verbose
                fprintf('      💾 Comprehensive gait plot saved: %s\n', plotPath);
            end
        end
        
    catch ME
        if verbose
            fprintf('      ❌ Error creating comprehensive gait plots: %s\n', ME.message);
        end
    end
end