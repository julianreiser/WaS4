function [heelStrikes, toeOffs] = detect_foot_events_fast(footZ, ankleAngle, timeVec, footSide, verbose)
% DETECT_FOOT_EVENTS_FAST - Ultra-sensitive detection for subtle motion capture signals
%
% Optimized for motion capture data with very small amplitude variations

heelStrikes = [];
toeOffs = [];

if length(footZ) < 50
    return;
end

% Estimate sampling rate
fs = 1 / mean(diff(timeVec));

if verbose
    fprintf('     🔍 %s foot: Signal range = %.4f units\n', footSide, range(footZ));
end

%% STEP 1: Enhanced preprocessing for subtle signals

% Remove overall trend (detrend)
footZ_detrended = detrend(footZ, 'linear');

% Apply bandpass filter to isolate gait frequencies (0.5-3 Hz for walking)
if fs > 6 % Need sufficient sampling rate
    [b, a] = butter(4, [0.5, 3] / (fs/2), 'bandpass');
    footZ_filtered = filtfilt(b, a, footZ_detrended);
else
    % Just smooth if sampling rate too low for bandpass
    footZ_filtered = movmean(footZ_detrended, 5);
end

% Calculate first and second derivatives for additional features
footZ_vel = gradient(footZ_filtered) * fs;
footZ_acc = gradient(footZ_vel) * fs;

if verbose
    fprintf('     📊 %s foot: Filtered signal range = %.4f units\n', footSide, range(footZ_filtered));
end

%% STEP 2: Adaptive threshold detection

% Use very sensitive thresholds based on signal statistics
signalStd = std(footZ_filtered);
signalRange = range(footZ_filtered);

% For heel strikes: use negative peaks (foot goes down)
% Use much smaller thresholds than before
minPeakHeight_HS = signalStd * 0.2; % Very sensitive
minPeakProminence_HS = signalStd * 0.1; % Very sensitive

% For toe offs: use positive peaks (foot goes up) 
minPeakHeight_TO = signalStd * 0.2;
minPeakProminence_TO = signalStd * 0.1;

% Minimum distance between steps (conservative estimate)
typical_step_time = 0.6; % seconds (100 steps/min)
minPeakDistance = round(typical_step_time * 0.7 * fs); % 70% of typical

if verbose
    fprintf('     🎯 %s foot: Thresholds - Height: %.4f, Prominence: %.4f, MinDist: %d samples\n', ...
            footSide, minPeakHeight_HS, minPeakProminence_HS, minPeakDistance);
end

%% STEP 3: Find peaks with very sensitive parameters

% Heel strikes: negative peaks in filtered signal
[~, hsLocs] = findpeaks(-footZ_filtered, ...
    'MinPeakHeight', minPeakHeight_HS, ...
    'MinPeakProminence', minPeakProminence_HS, ...
    'MinPeakDistance', minPeakDistance);

% Toe offs: positive peaks in filtered signal
[~, toLocs] = findpeaks(footZ_filtered, ...
    'MinPeakHeight', minPeakHeight_TO, ...
    'MinPeakProminence', minPeakProminence_TO, ...
    'MinPeakDistance', minPeakDistance);

if verbose
    fprintf('     🔍 %s foot: Initial detection - %d heel strikes, %d toe offs\n', ...
            footSide, length(hsLocs), length(toLocs));
end

%% STEP 4: Alternative method if peaks are still not found

if isempty(hsLocs) || isempty(toLocs)
    if verbose
        fprintf('     ⚠️ %s foot: Peak detection failed, trying alternative methods\n', footSide);
    end
    
    % Method A: Use zero-crossings of derivatives
    vel_zero_crossings = find(diff(sign(footZ_vel)) ~= 0);
    acc_zero_crossings = find(diff(sign(footZ_acc)) ~= 0);
    
    % Method B: Use local extrema without strict thresholds
    [hsLocs] = find_local_extrema(-footZ_filtered, minPeakDistance, 'min');
    [toLocs] = find_local_extrema(footZ_filtered, minPeakDistance, 'max');
    
    if verbose
        fprintf('     🔄 %s foot: Alternative detection - %d heel strikes, %d toe offs\n', ...
                footSide, length(hsLocs), length(toLocs));
    end
end

%% STEP 5: Temporal pattern validation

% Remove events that don't follow expected gait timing
validHS = [];
validTO = [];

% Expected gait cycle: HS -> TO (stance) -> HS (swing)
% Stance phase: typically 60% of gait cycle
% Swing phase: typically 40% of gait cycle

if ~isempty(hsLocs) && ~isempty(toLocs)
    
    % Start with heel strikes and find corresponding toe offs
    for i = 1:length(hsLocs)
        % Find toe offs after this heel strike
        laterTO = toLocs(toLocs > hsLocs(i));
        
        if ~isempty(laterTO)
            toIdx = laterTO(1);
            stanceTime = (toIdx - hsLocs(i)) / fs;
            
            % Very permissive stance time validation (0.1-1.2s)
            if stanceTime > 0.1 && stanceTime < 1.2
                validHS(end+1) = hsLocs(i);
                validTO(end+1) = toIdx;
            end
        end
    end
end

% If still no valid events, try the most permissive approach
if isempty(validHS)
    if verbose
        fprintf('     🆘 %s foot: Trying most permissive detection\n', footSide);
    end
    
    % Just find any local minima and maxima with minimal constraints
    validHS = find_local_extrema(-footZ_filtered, round(0.3*fs), 'min');
    validTO = find_local_extrema(footZ_filtered, round(0.3*fs), 'max');
    
    % Keep equal numbers
    minEvents = min(length(validHS), length(validTO));
    validHS = validHS(1:minEvents);
    validTO = validTO(1:minEvents);
end

%% STEP 6: Convert to time and final output

heelStrikes = timeVec(validHS);
toeOffs = timeVec(validTO);

if verbose
    fprintf('     ✅ %s foot: Final result - %d heel strikes, %d toe offs\n', ...
            footSide, length(heelStrikes), length(toeOffs));
    
    if ~isempty(heelStrikes)
        avgStepTime = mean(diff(heelStrikes));
        stepFreq = 1/avgStepTime;
        fprintf('     📈 %s foot: Avg step time = %.2fs, Freq = %.1f steps/min\n', ...
                footSide, avgStepTime, stepFreq*60);
    end
end

end

function extremaLocs = find_local_extrema(signal, minDist, type)
% FIND_LOCAL_EXTREMA - Find local minima or maxima without strict thresholds
%
% This is a fallback method for very subtle signals

extremaLocs = [];

if strcmp(type, 'min')
    % Find local minima
    for i = 2:(length(signal)-1)
        if signal(i) < signal(i-1) && signal(i) < signal(i+1)
            % Check if far enough from previous detection
            if isempty(extremaLocs) || (i - extremaLocs(end)) > minDist
                extremaLocs(end+1) = i;
            end
        end
    end
else % 'max'
    % Find local maxima
    for i = 2:(length(signal)-1)
        if signal(i) > signal(i-1) && signal(i) > signal(i+1)
            % Check if far enough from previous detection
            if isempty(extremaLocs) || (i - extremaLocs(end)) > minDist
                extremaLocs(end+1) = i;
            end
        end
    end
end

end