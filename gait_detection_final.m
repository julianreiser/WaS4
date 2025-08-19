function [heelStrikes, toeOffs, metrics] = gait_detection_final(EEG, foot, varargin)
% GAIT_DETECTION_FINAL - Final optimized gait event detection
%
% Detects heel strikes and toe-offs from XSens foot Z position data
% with parameters optimized for real walking data.
%
% USAGE: 
%   [heelStrikes, toeOffs, metrics] = gait_detection_final(EEG, 'left')
%   [heelStrikes, toeOffs, metrics] = gait_detection_final(EEG, 'right', 'Plot', true)
%   [heelStrikes, toeOffs, metrics] = gait_detection_final(EEG, 'left', 'TimeRange', [1400, 1600])
%
% INPUTS:
%   EEG      - EEGLAB structure with XSens foot data
%   foot     - 'left' or 'right'
%
% OUTPUTS:
%   heelStrikes - Array of heel strike times (seconds)
%   toeOffs     - Array of toe-off times (seconds)  
%   metrics     - Structure with gait timing metrics

% Parse inputs
p = inputParser;
addParameter(p, 'Plot', false, @islogical);
addParameter(p, 'TimeRange', [], @(x) length(x)==2);
addParameter(p, 'Verbose', true, @islogical);
parse(p, varargin{:});

plotResults = p.Results.Plot;
timeRange = p.Results.TimeRange;
verbose = p.Results.Verbose;

if verbose
    fprintf('\n=== GAIT EVENT DETECTION: %s FOOT ===\n', upper(foot));
end

% Find foot channel
channelLabels = {EEG.chanlocs.labels};
footPattern = sprintf('%s_Foot_Z', foot);
footIdx = find(strcmpi(channelLabels, footPattern), 1);

if isempty(footIdx)
    if verbose
        fprintf('❌ Channel %s not found!\n', footPattern);
        fprintf('Available foot channels: %s\n', ...
                strjoin(channelLabels(contains(channelLabels, 'Foot')), ', '));
    end
    heelStrikes = [];
    toeOffs = [];
    metrics = struct();
    return;
end

if verbose
    fprintf('✅ Channel found: %s\n', channelLabels{footIdx});
end

% Extract data
fs = EEG.srate;
timeVec = (0:EEG.pnts-1) / fs;
footZ = EEG.data(footIdx, :);

% Apply time range if specified
if ~isempty(timeRange)
    startIdx = round(timeRange(1) * fs) + 1;
    endIdx = round(timeRange(2) * fs);
    startIdx = max(1, startIdx);
    endIdx = min(length(footZ), endIdx);
    
    footZ = footZ(startIdx:endIdx);
    timeVec = timeVec(startIdx:endIdx);
    
    if verbose
        fprintf('🔍 Analyzing time range: %.1f - %.1f seconds (%d samples)\n', ...
                timeRange(1), timeRange(2), length(footZ));
    end
end

% Check data quality
validData = ~isnan(footZ) & ~isinf(footZ);
percentValid = 100 * sum(validData) / length(footZ);

if verbose
    fprintf('📊 Data quality: %.1f%% valid (%d/%d samples)\n', ...
            percentValid, sum(validData), length(footZ));
end

if sum(validData) < 100
    if verbose
        fprintf('❌ Insufficient valid data for analysis\n');
    end
    heelStrikes = [];
    toeOffs = [];
    metrics = struct();
    return;
end

% Handle NaN values by using only valid data segments
if percentValid < 90
    if verbose
        fprintf('⚠️  Using only valid data segments (%.1f%% of total)\n', percentValid);
    end
    validIdx = find(validData);
    footZ_work = footZ(validIdx);
    timeVec_work = timeVec(validIdx);
else
    footZ_work = footZ;
    timeVec_work = timeVec;
end

% Preprocessing: center and smooth lightly
footZ_work = footZ_work - mean(footZ_work, 'omitnan');
footZ_work = movmean(footZ_work, 3, 'omitnan'); % Very light smoothing

% Signal characteristics
signal_range = range(footZ_work);
signal_std = std(footZ_work, 'omitnan');

if verbose
    fprintf('📈 Signal characteristics:\n');
    fprintf('   Range: %.3f m, Std: %.3f m\n', signal_range, signal_std);
end

% Optimized parameters (based on successful detection)
prominence_factor = 0.01; % 1% of signal range
min_distance_sec = 0.3;   % 300ms minimum between events
min_distance = round(min_distance_sec * fs);

% Calculate thresholds
heel_strike_prominence = signal_range * prominence_factor * 0.3; % Lower for minima
toe_off_prominence = signal_range * prominence_factor; % Higher for maxima

if verbose
    fprintf('🎯 Detection parameters:\n');
    fprintf('   HS prominence: %.4f m (%.1f%% of range)\n', heel_strike_prominence, heel_strike_prominence/signal_range*100);
    fprintf('   TO prominence: %.4f m (%.1f%% of range)\n', toe_off_prominence, toe_off_prominence/signal_range*100);
    fprintf('   Min distance: %.1f seconds (%d samples)\n', min_distance_sec, min_distance);
end

% Detect heel strikes (foot minima)
[~, hsLocs] = findpeaks(-footZ_work, ...
    'MinPeakProminence', heel_strike_prominence, ...
    'MinPeakDistance', min_distance);

% Detect toe-offs (foot maxima)
[~, toLocs] = findpeaks(footZ_work, ...
    'MinPeakProminence', toe_off_prominence, ...
    'MinPeakDistance', min_distance);

% Convert to time
heelStrikes = timeVec_work(hsLocs);
toeOffs = timeVec_work(toLocs);

if verbose
    fprintf('\n📈 DETECTION RESULTS:\n');
    fprintf('   Heel Strikes: %d events\n', length(heelStrikes));
    fprintf('   Toe Offs: %d events\n', length(toeOffs));
    
    if ~isempty(heelStrikes) && ~isempty(toeOffs)
        all_events = [heelStrikes(:); toeOffs(:)]; % Ensure column vectors
        total_time = max(all_events) - min(all_events);
        step_rate = length(heelStrikes) / total_time * 60;
        fprintf('   Step rate: %.1f steps/minute\n', step_rate);
    end
end

% Calculate gait metrics
metrics = calculate_gait_metrics(heelStrikes, toeOffs, foot, verbose);

% Create visualization
if plotResults
    create_gait_validation_plot(footZ, timeVec, footZ_work, timeVec_work, ...
                               heelStrikes, toeOffs, foot, metrics);
end

end

function metrics = calculate_gait_metrics(heelStrikes, toeOffs, foot, verbose)
% Calculate comprehensive gait timing metrics

metrics = struct();
metrics.foot = foot;
metrics.n_heel_strikes = length(heelStrikes);
metrics.n_toe_offs = length(toeOffs);

% Step timing from heel strikes
if length(heelStrikes) > 1
    step_intervals = diff(heelStrikes);
    metrics.mean_step_time = mean(step_intervals);
    metrics.step_time_std = std(step_intervals);
    metrics.step_frequency = 1 / metrics.mean_step_time; % Hz
    metrics.cadence = metrics.step_frequency * 60; % steps/minute
    metrics.step_time_cov = metrics.step_time_std / metrics.mean_step_time * 100; % % coefficient of variation
end

% Stance and swing phase timing (if both events available)
if ~isempty(heelStrikes) && ~isempty(toeOffs)
    % Match heel strikes with subsequent toe-offs
    stance_times = [];
    swing_times = [];
    
    for i = 1:length(heelStrikes)
        % Find next toe-off after this heel strike
        next_toe_offs = toeOffs(toeOffs > heelStrikes(i));
        if ~isempty(next_toe_offs)
            stance_time = next_toe_offs(1) - heelStrikes(i);
            if stance_time > 0.1 && stance_time < 1.0 % Reasonable stance phase
                stance_times = [stance_times, stance_time];
            end
        end
        
        % Find next heel strike for swing time
        if i < length(heelStrikes)
            % Find toe-offs between this HS and next HS
            between_toe_offs = toeOffs(toeOffs > heelStrikes(i) & toeOffs < heelStrikes(i+1));
            if ~isempty(between_toe_offs)
                swing_time = heelStrikes(i+1) - between_toe_offs(end);
                if swing_time > 0.1 && swing_time < 0.8 % Reasonable swing phase
                    swing_times = [swing_times, swing_time];
                end
            end
        end
    end
    
    if ~isempty(stance_times)
        metrics.mean_stance_time = mean(stance_times);
        metrics.stance_time_std = std(stance_times);
        metrics.n_stance_phases = length(stance_times);
    end
    
    if ~isempty(swing_times)
        metrics.mean_swing_time = mean(swing_times);
        metrics.swing_time_std = std(swing_times);
        metrics.n_swing_phases = length(swing_times);
    end
    
    % Calculate stance/swing percentages
    if isfield(metrics, 'mean_stance_time') && isfield(metrics, 'mean_swing_time')
        total_cycle = metrics.mean_stance_time + metrics.mean_swing_time;
        metrics.stance_percentage = (metrics.mean_stance_time / total_cycle) * 100;
        metrics.swing_percentage = (metrics.mean_swing_time / total_cycle) * 100;
    end
end

if verbose && isfield(metrics, 'cadence')
    fprintf('\n📊 GAIT METRICS:\n');
    fprintf('   Cadence: %.1f steps/min\n', metrics.cadence);
    fprintf('   Step time: %.3f ± %.3f s (CV: %.1f%%)\n', ...
            metrics.mean_step_time, metrics.step_time_std, metrics.step_time_cov);
    
    if isfield(metrics, 'stance_percentage')
        fprintf('   Stance phase: %.1f%% (%.3f ± %.3f s)\n', ...
                metrics.stance_percentage, metrics.mean_stance_time, metrics.stance_time_std);
        fprintf('   Swing phase: %.1f%% (%.3f ± %.3f s)\n', ...
                metrics.swing_percentage, metrics.mean_swing_time, metrics.swing_time_std);
    end
end

end

function create_gait_validation_plot(footZ_orig, timeVec_orig, footZ_work, timeVec_work, ...
                                    heelStrikes, toeOffs, foot, metrics)
% Create comprehensive validation plot

figure('Position', [100, 100, 1400, 800], 'Name', sprintf('Gait Detection Results: %s Foot', upper(foot)));

% Plot 1: Original signal with events
subplot(2,1,1);
plot(timeVec_orig, footZ_orig, 'k-', 'LineWidth', 1, 'DisplayName', 'Original Signal');
hold on;

% Mark events on original signal
for i = 1:length(heelStrikes)
    hs_val = interp1(timeVec_orig, footZ_orig, heelStrikes(i), 'nearest', NaN);
    if ~isnan(hs_val)
        plot(heelStrikes(i), hs_val, 'rv', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
    end
end

for i = 1:length(toeOffs)
    to_val = interp1(timeVec_orig, footZ_orig, toeOffs(i), 'nearest', NaN);
    if ~isnan(to_val)
        plot(toeOffs(i), to_val, 'g^', 'MarkerSize', 8, 'MarkerFaceColor', 'green');
    end
end

xlabel('Time (seconds)');
ylabel('Foot Z Position (meters)');
title(sprintf('Original Signal: %s Foot with Detected Gait Events', upper(foot)));
legend('Original Signal', 'Heel Strikes', 'Toe Offs', 'Location', 'best');
grid on;

% Plot 2: Processed signal with events
subplot(2,1,2);
plot(timeVec_work, footZ_work, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Processed Signal');
hold on;

% Mark events on processed signal
for i = 1:length(heelStrikes)
    hs_val = interp1(timeVec_work, footZ_work, heelStrikes(i), 'nearest', NaN);
    if ~isnan(hs_val)
        plot(heelStrikes(i), hs_val, 'rv', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
    end
end

for i = 1:length(toeOffs)
    to_val = interp1(timeVec_work, footZ_work, toeOffs(i), 'nearest', NaN);
    if ~isnan(to_val)
        plot(toeOffs(i), to_val, 'g^', 'MarkerSize', 8, 'MarkerFaceColor', 'green');
    end
end

xlabel('Time (seconds)');
ylabel('Centered Foot Z Position (meters)');
title('Processed Signal with Gait Events');
legend('Processed Signal', 'Heel Strikes', 'Toe Offs', 'Location', 'best');
grid on;

% Add metrics text box
if isfield(metrics, 'cadence')
    metrics_text = {
        sprintf('GAIT METRICS (%s foot):', upper(foot)),
        sprintf('Events: %d HS, %d TO', metrics.n_heel_strikes, metrics.n_toe_offs),
        sprintf('Cadence: %.1f steps/min', metrics.cadence),
        sprintf('Step time: %.3f ± %.3f s', metrics.mean_step_time, metrics.step_time_std)
    };
    
    if isfield(metrics, 'stance_percentage')
        metrics_text{end+1} = sprintf('Stance: %.1f%%, Swing: %.1f%%', ...
                                     metrics.stance_percentage, metrics.swing_percentage);
    end
    
    annotation('textbox', [0.02, 0.02, 0.25, 0.2], 'String', metrics_text, ...
               'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', 'black');
end

end