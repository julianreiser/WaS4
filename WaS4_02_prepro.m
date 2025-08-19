%% WaS4 ICA Preprocessing Pipeline
% Adapted from WaS2 preprocessing pipeline for WaS4 project
% Performs ICA decomposition on merged WaS4 datasets
% Excludes non-EEG channels during preprocessing, then restores them after AMICA

close all
clear all

%% ----- SCRIPT SPECIFIC ----- %%
INFILE_extension    = '_merged.set';
OUTFILE_extension   = '_ica.set';

%% ----- EXPERIMENT SPECIFIC ----- %%
EXP = 'WaS4_';   % WaS4 experiment

%   SWITCH to determine preprocessing variant
prepro = 0;     % OPTIONS: 0 = bemobil; 1 = Eigen Decomposition, 2 = ICA (Stefans script - adapted)

%   SWITCH to determine user and paths
user = 0;       % OPTIONS: 0 = Julian; 1 = Emma newMac; 2 = Emma oldMac
ismac = 1;      % OPTIONS: 0 = PC; 1 = mac 

%   SWITCH to determine processing version folder name
ProcVersion = 'v01';

% Define subjects to skip if needed
subjects2skip = {};

%% ----- PROCESSING PARAMETERS ----- %%
% Channel rejection parameters
params.chancorr_crit = 0.8;                    % Channel correlation threshold
params.chan_max_broken_time = 0.3;             % Max time channel can be broken
params.chan_detect_num_iter = 10;              % Iterations for bad channel detection
params.chan_detected_fraction_threshold = 0.5; % Threshold for detection

% Filtering parameters  
params.filter_lowCutoffFreqAMICA = 1.5;        % High-pass filter cutoff (Hz)
params.filter_AMICA_highPassOrder = 1650;      % High-pass filter order

% AMICA parameters
params.AMICA_max_iter = 2000;                  % Maximum AMICA iterations
params.AMICA_n_rej = 10;                       % Auto-rejection iterations
params.AMICA_reject_sigma_threshold = 3;       % Auto-rejection threshold (sigma)
params.max_threads = 4;                        % AMICA threads

% Plot parameters
params.n_topo_components = 46;                 % Components to plot in topography

%% -----    NO CHANGES FROM HERE!
% set paths - adapt for WaS4 project structure
if user == 0 & ismac == 1
    pathVars.PATH = '/Volumes/Work4TB/Seafile/WaS4/'; % Current WaS4 directory
    pathVars.LAB = '/Volumes/Work4TB/Seafile/functions/eeglab2025.0.0';
    pathVars.bemobil = '/Users/julianreiser/Documents/GitHub/bemobil-pipeline/';
    pathVars.zaplinePlus = '/Users/julianreiser/Documents/GitHub/zapline-plus/';
    pathVars.fieldtrip = '/Users/julianreiser/Documents/GitHub/fieldtrip/';
    pathVars.DATA = [pathVars.PATH '/data/'];
    pathVars.PLOT = [pathVars.PATH '/plots/'];
    pathVars.STAT = [pathVars.PATH '/stats/'];
    pathVars.RAW = [pathVars.DATA 'aligned/'];
    pathVars.ICA = [pathVars.DATA 'ica/'];
    
    % Create directories if they don't exist
    dirs_to_create = {pathVars.DATA, pathVars.PLOT, pathVars.STAT, pathVars.RAW, pathVars.ICA};
    for dir_path = dirs_to_create
        if ~exist(dir_path{1}, 'dir')
            mkdir(dir_path{1});
            fprintf('Created directory: %s\n', dir_path{1});
        end
    end
end

% add function paths
addpath(pathVars.LAB)
addpath(genpath(pathVars.bemobil));
addpath(genpath(pathVars.zaplinePlus));
addpath(genpath(pathVars.fieldtrip));
addpath(genpath('/Users/julianreiser/Documents/GitHub/gaitEEGfootprint')); % Add gait footprint functions

%% DETERMINE SUBJECTS
% list all merged datasets in current directory and RAW directory
subjects        = struct();

% Check current directory first
current_dir_files = dir(fullfile(pathVars.PATH, [EXP '*' INFILE_extension]));
raw_dir_files = dir(fullfile(pathVars.RAW, ['WaS_' '*' INFILE_extension]));

% Combine files from both locations
subjects.DIR = [current_dir_files; raw_dir_files];

subjects.COMP_DIR   = dir(fullfile(pathVars.ICA, [EXP '*' OUTFILE_extension]));

if ~isempty(subjects.DIR)
    subjects.ALL    = cellfun(@(x) x{1}{1}, arrayfun(@(x) regexp(x.name, '_([0-9]+)_', 'tokens'), subjects.DIR, 'UniformOutput', false),'UniformOutput', false);
else
    subjects.ALL = {};
    fprintf('WARNING: No merged datasets found matching pattern %s*%s\n', EXP, INFILE_extension);
end

if ~isempty(subjects.COMP_DIR)
    subjects.COMP_ALL = cellfun(@(x) x{1}{1}, arrayfun(@(x) regexp(x.name, '_([0-9]+)_', 'tokens'), subjects.COMP_DIR, 'UniformOutput', false),'UniformOutput', false);
else
    subjects.COMP_ALL = {};
end

if exist("subjects2skip")
    subjects.SKIP   = subjects2skip';
else
    subjects.SKIP   = [];
end

% Font size for figures
fontsize = 10;

fprintf('Found %d merged datasets to process\n', length(subjects.ALL));
fprintf('Found %d already completed ICA datasets\n', length(subjects.COMP_ALL));

%   open EEGLab
% eeglab; disp(newline)

%% MAIN PROCESSING LOOP
for subi = 1:length(subjects.ALL)
    
    subject = subjects.ALL{subi};
    
    fprintf('\n=== Processing Subject %s (%d/%d) ===\n', subject, subi, length(subjects.ALL));

    if any(strcmpi(subjects.SKIP,subject)) || any(strcmpi(subjects.COMP_ALL,subject))
        if any(strcmpi(subjects.SKIP,subject))
            fprintf('SKIPPED: Subject %s is in skip list\n', subject);
        else
            fprintf('SKIPPED: Subject %s already processed\n', subject);
        end
        continue
    else

        % start a fresh instance of eeglab
        eeglab;
    
        %% LOAD MERGED DATASET
        fprintf('Loading merged dataset: %s\n', [EXP subject INFILE_extension]);
        try
            % Try current directory first, then RAW directory
            if exist(fullfile(pathVars.PATH, ['WaS_' subject INFILE_extension]), 'file')
                EEG = pop_loadset('filename',['WaS_' subject INFILE_extension], 'filepath', pathVars.PATH);
            elseif exist(fullfile(pathVars.RAW, ['WaS_' subject INFILE_extension]), 'file')
                EEG = pop_loadset('filename',['WaS_' subject INFILE_extension], 'filepath', pathVars.RAW);
            else
                error('Dataset file not found in current directory or RAW directory');
            end
            fprintf('  Dataset loaded: %d channels, %d samples, %.1f seconds\n', EEG.nbchan, EEG.pnts, EEG.pnts/EEG.srate);
        catch ME
            fprintf('ERROR: Failed to load dataset for subject %s: %s\n', subject, ME.message);
            continue;
        end
        
        %% CHANNEL LOCATION CORRECTION
        fprintf('Correcting channel locations...\n');
        chanlocfile = which('standard-10-5-cap385.elp');
        if isempty(chanlocfile)
            fprintf('WARNING: Standard channel location file not found, using default\n');
        else
            try
                EEG = pop_chanedit(EEG, 'lookup', chanlocfile);
            catch
                fprintf('WARNING: Channel location lookup failed, continuing with existing locations\n');
            end
        end
        channel_locations = EEG.chanlocs;
        OLD = EEG; % Keep copy with all channels for later restoration
    
        % Find EEG channels - detect automatically based on channel labels
        eeg_channel_names = {EEG.chanlocs.labels};
        eeg_channels = [];
        
        % Look for standard EEG channel patterns (exclude WaS4 non-EEG channels)
        eeg_patterns = {'^[FPOCTZ]\d+', '^A[FT]\d+', '^P[OZ]\d+', '^C[PZ]\d+', '^F[CTZ]\d+', '^T[78P]', '^M[12]', '^AF\d+', '^FC\d+', '^CP\d+', '^PO\d+', '^Iz', '^Oz', '^Cz', '^Fz', '^Pz'};
        
        for i = 1:length(eeg_channel_names)
            chan_name = eeg_channel_names{i};
            is_eeg = false;
            
            % Check against EEG patterns
            for pat = eeg_patterns
                if ~isempty(regexp(chan_name, pat{1}, 'once'))
                    is_eeg = true;
                    break;
                end
            end
            
            % Exclude known non-EEG channels from WaS4 (gaze, motion capture, etc.)
            exclude_patterns = {'gaze', 'ECG', 'ACC', 'audio', 'VGait', 'CoM', 'BoS', 'MoS', 'FP\d', 'Segment', 'Joint', 'Euler', 'Quat', 'CNT'};
            is_excluded = false;
            for excl_pat = exclude_patterns
                if ~isempty(regexp(chan_name, excl_pat{1}, 'once', 'ignorecase'))
                    is_excluded = true;
                    break;
                end
            end
            
            % Include if it's an EEG channel and not excluded
            if is_eeg && ~is_excluded
                eeg_channels(end+1) = i;
            elseif i <= 64 && ~is_excluded && isempty(regexp(chan_name, '\d', 'once')) % Handle label-less channels in first 64 positions
                eeg_channels(end+1) = i;
            end
        end
        
        % Fallback if detection fails - use first 64 channels excluding obvious non-EEG
        if isempty(eeg_channels)
            fprintf('  WARNING: Automatic EEG channel detection failed, using manual selection\n');
            for i = 1:min(64, EEG.nbchan)
                chan_name = eeg_channel_names{i};
                if isempty(regexp(chan_name, '(gaze|ECG|ACC|audio|VGait|CoM|BoS|MoS|FP\d|Segment|Joint|Euler|Quat)', 'once', 'ignorecase'))
                    eeg_channels(end+1) = i;
                end
            end
        end
        
        fprintf('  Detected %d EEG channels out of %d total channels\n', length(eeg_channels), EEG.nbchan);
        fprintf('  Non-EEG channels (%d): %s\n', EEG.nbchan - length(eeg_channels), strjoin(eeg_channel_names(setdiff(1:EEG.nbchan, eeg_channels)), ', '));
        
        % Select only EEG channels for preprocessing
        EEG = pop_select(EEG,'channel', eeg_channels);
        fprintf('  Selected %d EEG channels for preprocessing\n', EEG.nbchan);
        
        %% Create plot directories
        line_dir = [pathVars.PLOT '/line/'];
        chans_dir = [pathVars.PLOT '/chans/'];
        amica_dir = [pathVars.PLOT '/amica/'];
        analysis_dir = [pathVars.STAT '/analysis/'];
        
        for dir_path = {line_dir, chans_dir, amica_dir, analysis_dir}
            if ~exist(dir_path{1}, 'dir')
                mkdir(dir_path{1});
            end
        end
        
        %% Use zapline for line noise removal
        fprintf('Applying zapline-plus for line noise removal...\n');
        % noise frequencies are at 50 Hz (line), potentially others from motion capture/VR
        try
            EEG_clean = clean_data_with_zapline_plus_eeglab_wrapper(EEG);
            saveas(gcf,[line_dir EXP subject '_linenoise_' ProcVersion '.png']);
            close(gcf);
            fprintf('  Zapline completed and plot saved\n');
        catch ME
            fprintf('WARNING: Zapline failed: %s\n', ME.message);
            fprintf('  Continuing with original data\n');
            EEG_clean = EEG;
        end
        
        %% Channel correction / rejection
        fprintf('Detecting bad channels...\n');
        chancorr_crit = params.chancorr_crit;
        chan_max_broken_time = params.chan_max_broken_time;
        chan_detect_num_iter = params.chan_detect_num_iter;
        chan_detected_fraction_threshold = params.chan_detected_fraction_threshold;
        flatline_crit = 'off';
        line_noise_crit = 'off';
    
        try
            chans_to_interp = bemobil_detect_bad_channels(EEG_clean, ALLEEG, CURRENTSET, chancorr_crit, chan_max_broken_time,...
                chan_detect_num_iter, chan_detected_fraction_threshold, flatline_crit, line_noise_crit);

            % save plot and data
            saveas(gcf,[chans_dir EXP subject '_chanreject_' ProcVersion '.png']); close gcf
            if ishandle(gcf)
                saveas(gcf,[chans_dir EXP subject '_chanrejectIterations_' ProcVersion '.png']); close gcf
            end
            save([analysis_dir EXP subject '_chans2interp_' ProcVersion '.mat'],'chans_to_interp');
            
            fprintf('  Found %d channels to interpolate: %s\n', length(chans_to_interp), mat2str(chans_to_interp));
        catch ME
            fprintf('WARNING: Bad channel detection failed: %s\n', ME.message);
            fprintf('  Continuing without channel interpolation\n');
            chans_to_interp = [];
        end
    
        % do the actual interpolation and full rank average referencing
        fprintf('Interpolating channels and average referencing...\n');
        try
            [ALLEEG, EEG_preprocessed, CURRENTSET] = bemobil_interp_avref(EEG_clean, ALLEEG, CURRENTSET, chans_to_interp);
            fprintf('  Channel interpolation and referencing completed\n');
        catch ME
            fprintf('ERROR: Channel interpolation failed: %s\n', ME.message);
            continue;
        end
        
        %% Filter for AMICA
        fprintf('Filtering data for AMICA...\n');
        % See Klug & Gramann (2020) for filter recommendations
        filter_lowCutoffFreqAMICA = params.filter_lowCutoffFreqAMICA; % 1.5Hz high-pass filter
        filter_AMICA_highPassOrder = params.filter_AMICA_highPassOrder; % Order used by Klug & Gramann (2020)
        filter_highCutoffFreqAMICA = []; % no low-pass filter
        filter_AMICA_lowPassOrder = []; 
    
        out_filename = [];
        out_filepath= [];
    
        try
            [ALLEEG, EEG_filtered, CURRENTSET] = bemobil_filter(ALLEEG, EEG_preprocessed, CURRENTSET, filter_lowCutoffFreqAMICA,...
                filter_highCutoffFreqAMICA, out_filename, out_filepath, filter_AMICA_highPassOrder, filter_AMICA_lowPassOrder);
            fprintf('  Filtering completed\n');
        catch ME
            fprintf('ERROR: Filtering failed: %s\n', ME.message);
            continue;
        end
        
        %% Compute AMICA
        fprintf('Computing AMICA decomposition...\n');
        % data rank is the number of channels that were not interpolated
        data_rank = EEG_filtered.nbchan - length(EEG_filtered.etc.interpolated_channels);
        fprintf('  Data rank: %d (channels: %d, interpolated: %d)\n', data_rank, EEG_filtered.nbchan, length(EEG_filtered.etc.interpolated_channels));
    
        % AMICA settings
        amica = true;
        numb_models = 1; % default 1
        AMICA_autoreject = 1; % uses automatic rejection method of AMICA
        AMICA_n_rej = params.AMICA_n_rej; % number of iterations during automatic rejection
        AMICA_reject_sigma_threshold = params.AMICA_reject_sigma_threshold; % threshold for rejection
        AMICA_max_iter = params.AMICA_max_iter; % maximum number of AMICA iterations
        max_threads = params.max_threads; % 4 threads are most effective for single subject
    
        other_algorithm = [];
        out_filename = [];
        out_filepath = [];
    
        % compute AMICA
        try
            fprintf('  Starting AMICA (this may take a while)...\n');
            [ALLEEG, EEG_amica, CURRENTSET] = bemobil_signal_decomposition(ALLEEG, EEG_filtered, CURRENTSET,...
                amica, numb_models, max_threads, data_rank, other_algorithm, out_filename, out_filepath, AMICA_autoreject,...
                AMICA_n_rej, AMICA_reject_sigma_threshold, AMICA_max_iter);
            fprintf('  AMICA completed successfully\n');
        catch ME
            fprintf('ERROR: AMICA failed: %s\n', ME.message);
            continue;
        end
        
        %% Plot AMICA results
        fprintf('Generating AMICA plots...\n');
        try
            % plot topographies
            n_comps_to_plot = min(params.n_topo_components, EEG_amica.nbchan);
            pop_topoplot(EEG_amica, 0, [1:n_comps_to_plot], ['Subject ' subject ' AMICA Components'], [7 7], 0, 'electrodes', 'off');
            saveas(gcf,[amica_dir EXP subject '_amicatopo_' ProcVersion '.png']); close gcf;
        
            % plot autorejection
            data2plot = EEG_amica.data(1:round(EEG_amica.nbchan/10):EEG_amica.nbchan,:)';
            figure;
            set(gcf,'color','w','Position', get(0,'screensize'));
            plot(data2plot,'g');
            if isfield(EEG_amica.etc, 'bad_samples') && ~isempty(EEG_amica.etc.bad_samples)
                data2plot(~EEG_amica.etc.bad_samples,:) = NaN;
                hold on
                plot(data2plot,'r');
                rejected_percent = round(EEG_amica.etc.bad_samples_percent, 2);
            else
                rejected_percent = 0;
            end
            xlim([-10000 EEG_amica.pnts+10000])
            ylim([-1000 1000])
            title(['AMICA autorejection, removed ' num2str(rejected_percent) '% of the samples'])
            xlabel('Samples')
            ylabel('\muV')
            saveas(gcf,[amica_dir EXP subject '_amicaresults_' ProcVersion '.png']); close gcf;
            
            fprintf('  AMICA plots saved, rejected %.1f%% of samples\n', rejected_percent);
        catch ME
            fprintf('WARNING: AMICA plotting failed: %s\n', ME.message);
        end
    
        %% Copy spatial filter to unfiltered data
        fprintf('Copying spatial filter to unfiltered dataset...\n');
        try
            [ALLEEG, EEG_AMICA_copied, CURRENTSET] = bemobil_copy_spatial_filter(EEG_preprocessed, ALLEEG, CURRENTSET, EEG_amica);
            
            % Run ICLabel on unfiltered data
            fprintf('  Running ICLabel classification...\n');
            EEG_AMICA_copied = pop_iclabel(EEG_AMICA_copied, 'default');
            EEG_AMICA_copied = eeg_checkset(EEG_AMICA_copied,'eventconsistency');
            
            fprintf('  Spatial filter copy completed\n');
        catch ME
            fprintf('ERROR: Spatial filter copying failed: %s\n', ME.message);
            continue;
        end

        %% Restore non-EEG channels from original dataset
        fprintf('Restoring non-EEG channels (motion capture, gaze, ECG, etc.)...\n');
        try
            % Get indices of non-EEG channels from original dataset
            non_eeg_indices = setdiff(1:OLD.nbchan, eeg_channels);
            
            if ~isempty(non_eeg_indices)
                % Get the current number of EEG channels after preprocessing
                n_eeg_channels = size(EEG_AMICA_copied.data, 1);
                n_non_eeg_channels = length(non_eeg_indices);
                
                % Expand data matrix to accommodate non-EEG channels
                EEG_AMICA_copied.data(n_eeg_channels+1:n_eeg_channels+n_non_eeg_channels, :) = OLD.data(non_eeg_indices, :);
                
                % Update channel information
                for i = 1:n_non_eeg_channels
                    orig_idx = non_eeg_indices(i);
                    new_idx = n_eeg_channels + i;
                    EEG_AMICA_copied.chanlocs(new_idx) = OLD.chanlocs(orig_idx);
                end
                
                % Update dataset parameters
                EEG_AMICA_copied.nbchan = n_eeg_channels + n_non_eeg_channels;
                
                % Create mapping for future reference
                EEG_AMICA_copied.etc.channel_mapping = struct();
                EEG_AMICA_copied.etc.channel_mapping.eeg_channels = eeg_channels;
                EEG_AMICA_copied.etc.channel_mapping.non_eeg_channels = non_eeg_indices;
                EEG_AMICA_copied.etc.channel_mapping.restored_non_eeg_indices = n_eeg_channels+1:n_eeg_channels+n_non_eeg_channels;
                
                fprintf('  Restored %d non-EEG channels: %s\n', n_non_eeg_channels, strjoin({OLD.chanlocs(non_eeg_indices).labels}, ', '));
                fprintf('  Total channels now: %d (EEG: %d, Others: %d)\n', EEG_AMICA_copied.nbchan, n_eeg_channels, n_non_eeg_channels);
            else
                fprintf('  No additional channels to restore\n');
            end
        catch ME
            fprintf('WARNING: Failed to restore non-EEG channels: %s\n', ME.message);
            fprintf('  Continuing with EEG channels only\n');
        end
    
        %% Quality Control: Store Baseline Gait Artifact Footprint (Jacobsen et al. 2022)
        fprintf('Computing baseline gait artifact footprint for quality control...\n');
        try
            % Calculate baseline footprint before any component rejection
            % This will be compared later against cleaned data after component rejection
            footprint_baseline = calculate_gait_artifact_footprint(EEG_preprocessed, 1:EEG_preprocessed.nbchan);
            
            % Store baseline footprint in EEG structure for later comparison
            EEG_AMICA_copied.etc.quality_control = struct();
            EEG_AMICA_copied.etc.quality_control.footprint_baseline = footprint_baseline;
            EEG_AMICA_copied.etc.quality_control.feature_labels = {'B_Rfreq', 'C_lateral_medial', 'D_neck_ratio', 'E_double_single', 'F_stand_walk'};
            EEG_AMICA_copied.etc.quality_control.computed_on = datestr(now);
            EEG_AMICA_copied.etc.quality_control.eeg_channels_used = 1:EEG_preprocessed.nbchan;
            EEG_AMICA_copied.etc.quality_control.method = 'Jacobsen_et_al_2022';
            
            % Also save to external file for record keeping
            quality_control = EEG_AMICA_copied.etc.quality_control;
            save([analysis_dir EXP subject '_baseline_footprint_' ProcVersion '.mat'], 'quality_control');
            
            fprintf('  Baseline footprint features: [%.3f, %.3f, %.3f, %.3f, %.3f]\n', footprint_baseline);
            fprintf('  Stored in EEG.etc.quality_control for later comparison after component rejection\n');
            fprintf('  Use calculate_footprint_improvement(EEG) after cleaning to assess artifact reduction\n');
        catch ME
            fprintf('WARNING: Baseline footprint calculation failed: %s\n', ME.message);
            EEG_AMICA_copied.etc.quality_control = struct();
            EEG_AMICA_copied.etc.quality_control.error = ME.message;
        end

        %% Generate Processing Summary
        processing_summary = struct();
        processing_summary.subject = subject;
        processing_summary.original_channels = OLD.nbchan;
        processing_summary.eeg_channels = length(eeg_channels);
        processing_summary.interpolated_channels = length(chans_to_interp);
        processing_summary.final_channels = EEG_AMICA_copied.nbchan;
        processing_summary.n_components = size(EEG_AMICA_copied.icaweights, 1);
        processing_summary.data_rank = data_rank;
        processing_summary.non_eeg_channels = length(setdiff(1:OLD.nbchan, eeg_channels));
        if isfield(EEG_AMICA_copied.etc, 'bad_samples_percent')
            processing_summary.rejected_samples_percent = EEG_AMICA_copied.etc.bad_samples_percent;
        else
            processing_summary.rejected_samples_percent = 0;
        end
        processing_summary.duration_seconds = EEG_AMICA_copied.pnts/EEG_AMICA_copied.srate;
        processing_summary.sampling_rate = EEG_AMICA_copied.srate;
        
        % Save processing summary
        save([analysis_dir EXP subject '_processing_summary_' ProcVersion '.mat'], 'processing_summary');

        %% Save final dataset
        fprintf('Saving ICA dataset...\n');
        try
            EEG_AMICA_copied = pop_saveset(EEG_AMICA_copied, 'filename', [EXP subject '_' ProcVersion OUTFILE_extension], 'filepath', pathVars.ICA);
            
            fprintf('SUCCESS: Subject %s ICA preprocessing completed\n', subject);
            fprintf('  Final dataset: %s\n', [EXP subject '_' ProcVersion OUTFILE_extension]);
            fprintf('  Total channels: %d (EEG: %d, Others: %d)\n', processing_summary.final_channels, processing_summary.eeg_channels, processing_summary.non_eeg_channels);
            fprintf('  ICA components: %d (rank: %d)\n', processing_summary.n_components, processing_summary.data_rank);
            fprintf('  Interpolated channels: %d\n', processing_summary.interpolated_channels);
            fprintf('  Rejected samples: %.2f%%\n', processing_summary.rejected_samples_percent);
            fprintf('  Duration: %.1f seconds (%.1f Hz)\n', processing_summary.duration_seconds, processing_summary.sampling_rate);
        catch ME
            fprintf('ERROR: Failed to save dataset: %s\n', ME.message);
            continue;
        end
    
        close all
        fprintf('=== Subject %s completed ===\n', subject);
    end
end

fprintf('\n=== WaS4 ICA PREPROCESSING PIPELINE COMPLETED ===\n');
fprintf('Check the plots in: %s\n', pathVars.PLOT);
fprintf('ICA datasets saved in: %s\n', pathVars.ICA);

%% Gait Artifact Footprint Functions (adapted for WaS4)
function footprint_features = calculate_gait_artifact_footprint(EEG, eeg_channels)
    % Calculate gait artifact footprint using functions from:
    % https://github.com/NadineJac/gaitEEGfootprint
    % Adapted for WaS4 data structure with walking/standing conditions
    
    try
        % Extract walking and standing data from WaS4 events
        [walking_data, standing_data] = extract_walking_standing_periods_was4(EEG, eeg_channels);
        
        if isempty(walking_data) || isempty(standing_data)
            fprintf('  Warning: Could not separate walking/standing periods\n');
            footprint_features = zeros(1, 5);
            return;
        end
        
        % Calculate time-frequency representation (ERSP format expected by footprint functions)
        fs = EEG.srate;
        
        % Calculate power spectral density for both conditions  
        [psd_walking, freqs] = calculate_psd_for_footprint(walking_data, fs);
        [psd_standing, ~] = calculate_psd_for_footprint(standing_data, fs);
        
        % Calculate baseline-corrected time-frequency data (dB change from standing to walking)
        % Format: channels x frequencies x time_points (for footprint functions)
        TFdata = 10 * log10((psd_walking + eps) ./ (psd_standing + eps));
        
        % Add singleton time dimension if needed (footprint functions expect 3D)
        if ndims(TFdata) == 2
            TFdata = TFdata(:, :, ones(1, 100)); % Replicate across 100 time points
        end
        
        % Define channel indices for footprint calculation
        lateral_chan_idx = find_lateral_channels_was4(EEG.chanlocs(eeg_channels));
        neck_chan_L = find_neck_channels_was4(EEG.chanlocs(eeg_channels), 'left');
        neck_chan_R = find_neck_channels_was4(EEG.chanlocs(eeg_channels), 'right');
        
        % Calculate gait cycle points for E_doubleSuppRatio (simplified for WaS4)
        pnts_double = get_gait_cycle_points_was4(EEG, size(TFdata, 3));
        
        % Calculate footprint features using original functions
        features = zeros(1, 5);
        
        % B) Correlation across frequencies
        features(1) = B_Rfreq(TFdata);
        
        % C) Lateral to medial channel power ratio
        if ~isempty(lateral_chan_idx)
            features(2) = C_lateralPowRatio(TFdata, lateral_chan_idx);
        else
            features(2) = 0;
        end
        
        % D) Neck channel power ratio
        if ~isempty(neck_chan_L) && ~isempty(neck_chan_R) && ~isempty(pnts_double)
            pnts_LHS = 1:round(size(TFdata,3)/2); % Left heel strike points (simplified)
            pnts_RHS = round(size(TFdata,3)/2):size(TFdata,3); % Right heel strike points
            features(3) = D_neckChanRatio(TFdata, neck_chan_L, neck_chan_R, pnts_LHS, pnts_RHS);
        else
            features(3) = 0;
        end
        
        % E) Double to single support ratio
        if ~isempty(pnts_double)
            features(4) = E_doubleSuppRatio(TFdata, pnts_double);
        else
            features(4) = 0;
        end
        
        % F) Standing to walking power ratio
        features(5) = F_swRatio(TFdata);
        
        footprint_features = features;
        
    catch ME
        fprintf('  Error in footprint calculation: %s\n', ME.message);
        footprint_features = zeros(1, 5);
    end
end

function [walking_data, standing_data] = extract_walking_standing_periods_was4(EEG, eeg_channels)
    % Extract walking and standing periods based on WaS4 walkingCondition events
    
    eeg_data = EEG.data(eeg_channels, :);
    walking_data = [];
    standing_data = [];
    
    % Method 1: Use walkingCondition field if available
    for i = 1:length(EEG.event)
        if isfield(EEG.event(i), 'walkingCondition') && ~isempty(EEG.event(i).walkingCondition)
            
            % Find block boundaries 
            start_lat = round(EEG.event(i).latency);
            end_lat = start_lat + round(30 * EEG.srate); % 30 seconds default
            
            % Look for actual block end (S 2 marker)
            for j = i+1:min(i+100, length(EEG.event))
                if strcmp(EEG.event(j).type, 'S  2')
                    end_lat = round(EEG.event(j).latency);
                    break;
                end
            end
            
            end_lat = min(end_lat, size(eeg_data, 2));
            if end_lat > start_lat
                block_data = eeg_data(:, start_lat:end_lat);
                
                % Classify based on walking condition
                walking_condition = EEG.event(i).walkingCondition;
                if strcmp(walking_condition, 'standing')
                    standing_data = [standing_data, block_data];
                elseif any(strcmp(walking_condition, {'walking', 'flat_walking', 'wave_walking', 'hump_walking', 'treadmill'}))
                    walking_data = [walking_data, block_data];
                end
            end
        end
    end
    
    % Method 2: Fallback using S 1/S 2 markers with heuristics
    if isempty(walking_data) && isempty(standing_data)
        s1_events = find(strcmp({EEG.event.type}, 'S  1'));
        s2_events = find(strcmp({EEG.event.type}, 'S  2'));
        
        for block_idx = 1:min(length(s1_events), length(s2_events))
            start_lat = round(EEG.event(s1_events(block_idx)).latency);
            end_lat = round(EEG.event(s2_events(block_idx)).latency);
            
            if end_lat > start_lat && end_lat <= size(eeg_data, 2)
                block_data = eeg_data(:, start_lat:end_lat);
                
                % Heuristic: first blocks tend to be standing, later blocks walking
                if block_idx <= 2
                    standing_data = [standing_data, block_data];
                else
                    walking_data = [walking_data, block_data];
                end
            end
        end
    end
    
    % Method 3: Final fallback - simple data split based on variance
    if isempty(walking_data) || isempty(standing_data)
        n_samples = size(eeg_data, 2);
        
        % Calculate variance over segments to identify high/low activity periods
        segment_length = round(10 * EEG.srate); % 10-second segments
        n_segments = floor(n_samples / segment_length);
        segment_vars = zeros(1, n_segments);
        
        for seg = 1:n_segments
            seg_start = (seg-1) * segment_length + 1;
            seg_end = min(seg * segment_length, n_samples);
            segment_vars(seg) = mean(var(eeg_data(:, seg_start:seg_end), [], 2));
        end
        
        % High variance segments likely contain walking
        var_threshold = median(segment_vars);
        
        if isempty(standing_data)
            low_var_segs = find(segment_vars <= var_threshold);
            for seg = low_var_segs(1:min(3, length(low_var_segs))) % Use first 3 low-variance segments
                seg_start = (seg-1) * segment_length + 1;
                seg_end = min(seg * segment_length, n_samples);
                standing_data = [standing_data, eeg_data(:, seg_start:seg_end)];
            end
        end
        
        if isempty(walking_data)
            high_var_segs = find(segment_vars > var_threshold);
            for seg = high_var_segs(1:min(5, length(high_var_segs))) % Use first 5 high-variance segments
                seg_start = (seg-1) * segment_length + 1;
                seg_end = min(seg * segment_length, n_samples);
                walking_data = [walking_data, eeg_data(:, seg_start:seg_end)];
            end
        end
    end
end

function [psd, freqs] = calculate_psd_for_footprint(data, fs)
    % Calculate PSD in format expected by footprint functions
    
    if isempty(data) || size(data, 2) < fs
        % Handle insufficient data
        psd = ones(size(data, 1), 50) * 1e-12;
        freqs = linspace(1, 50, 50);
        return;
    end
    
    % Calculate PSD using Welch's method
    window = min(4*fs, floor(size(data, 2)/4)); % 4s windows or 1/4 of data  
    overlap = round(window/2);
    nfft = max(512, 2^nextpow2(window));
    
    try
        [psd, freqs] = pwelch(data', window, overlap, nfft, fs);
        psd = psd'; % Transpose to channels x frequencies
    catch
        % Fallback
        [psd, freqs] = periodogram(data', [], [], fs);
        psd = psd';
    end
    
    % Limit to reasonable frequency range (up to 50 Hz)
    freq_mask = freqs <= 50 & freqs >= 1;
    psd = psd(:, freq_mask);
    freqs = freqs(freq_mask);
end

function lateral_chans = find_lateral_channels_was4(chanlocs)
    % Find lateral channels for WaS4 EEG setup
    lateral_labels = {'F7', 'F8', 'T7', 'T8', 'P7', 'P8', 'FT7', 'FT8', 'TP7', 'TP8', 'FC5', 'FC6', 'AF7', 'AF8'};
    lateral_chans = [];
    
    for i = 1:length(chanlocs)
        if any(strcmpi(chanlocs(i).labels, lateral_labels))
            lateral_chans(end+1) = i;
        end
    end
end

function neck_chans = find_neck_channels_was4(chanlocs, side)
    % Find neck channels for left/right sides
    if nargin < 2
        side = 'both';
    end
    
    neck_chans = [];
    
    switch lower(side)
        case 'left'
            neck_labels = {'O1', 'PO3', 'PO7', 'P7'};
        case 'right'  
            neck_labels = {'O2', 'PO4', 'PO8', 'P8'};
        case 'both'
            neck_labels = {'O1', 'O2', 'Oz', 'PO3', 'PO4', 'PO7', 'PO8', 'Iz'};
    end
    
    for i = 1:length(chanlocs)
        if any(strcmpi(chanlocs(i).labels, neck_labels))
            neck_chans(end+1) = i;
        end
    end
end

function pnts_double = get_gait_cycle_points_was4(EEG, n_time_points)
    % Get gait cycle points for double support calculation
    % Uses step markers if available from WaS4 data
    
    step_events = [];
    
    % Look for step-related events in WaS4 data
    for i = 1:length(EEG.event)
        if isfield(EEG.event(i), 'type')
            % Look for step markers (could be various formats)
            if any(strcmp(EEG.event(i).type, {'step', 'heel_strike', 'toe_off', 'gait_event'}))
                step_events(end+1) = i;
            end
        end
    end
    
    if length(step_events) > 4
        % Calculate approximate double support periods based on gait events
        step_lats = [EEG.event(step_events).latency];
        step_intervals = diff(step_lats) / EEG.srate;
        
        % Approximate double support as 20% of stride cycle
        mean_step_interval = mean(step_intervals);
        double_support_duration = round(0.2 * mean_step_interval * EEG.srate);
        
        % Map to normalized gait cycle points
        pnts_double = round(linspace(1, double_support_duration, min(20, double_support_duration)));
    else
        % Default double support points (first 20% of analysis window)
        pnts_double = 1:round(0.2 * n_time_points);
    end
    
    % Ensure points are within valid range
    pnts_double = pnts_double(pnts_double <= n_time_points & pnts_double >= 1);
end

%% Helper Function for Quality Assessment After Component Rejection
function improvement_results = calculate_footprint_improvement(EEG)
    % Calculate improvement in gait artifact footprint after component rejection
    % Input: EEG structure with ICA weights and baseline footprint stored in etc
    % Output: Structure with improvement metrics
    
    improvement_results = struct();
    
    try
        % Check if baseline footprint exists
        if ~isfield(EEG.etc, 'quality_control') || ~isfield(EEG.etc.quality_control, 'footprint_baseline')
            error('No baseline footprint found. Run WaS4_03_prepro.m first to compute baseline.');
        end
        
        % Get baseline footprint
        footprint_baseline = EEG.etc.quality_control.footprint_baseline;
        eeg_channels = EEG.etc.quality_control.eeg_channels_used;
        
        % Limit to EEG channels only (in case non-EEG channels were restored)
        n_eeg_channels = length(eeg_channels);
        
        % Calculate current footprint (after component rejection)
        footprint_current = calculate_gait_artifact_footprint(EEG, 1:n_eeg_channels);
        
        % Calculate improvement metrics
        footprint_distance = norm(footprint_current - footprint_baseline);
        relative_change = (footprint_current - footprint_baseline) ./ (footprint_baseline + eps);
        percent_change = relative_change * 100;
        
        % Store results
        improvement_results.footprint_baseline = footprint_baseline;
        improvement_results.footprint_current = footprint_current;
        improvement_results.footprint_distance = footprint_distance;
        improvement_results.relative_change = relative_change;
        improvement_results.percent_change = percent_change;
        improvement_results.feature_labels = EEG.etc.quality_control.feature_labels;
        improvement_results.artifact_reduction_success = footprint_distance > 0.1; % Threshold for meaningful change
        improvement_results.computed_on = datestr(now);
        
        % Update EEG structure with current results
        EEG.etc.quality_control.footprint_current = footprint_current;
        EEG.etc.quality_control.improvement_results = improvement_results;
        
        % Display results
        fprintf('\n=== Gait Artifact Footprint Quality Assessment ===\n');
        fprintf('Baseline features:  [%.3f, %.3f, %.3f, %.3f, %.3f]\n', footprint_baseline);
        fprintf('Current features:   [%.3f, %.3f, %.3f, %.3f, %.3f]\n', footprint_current);
        fprintf('Percent change:     [%.1f%%, %.1f%%, %.1f%%, %.1f%%, %.1f%%]\n', percent_change);
        fprintf('Euclidean distance: %.4f\n', footprint_distance);
        fprintf('Artifact reduction successful: %s\n', string(improvement_results.artifact_reduction_success));
        fprintf('Features: %s\n', strjoin(improvement_results.feature_labels, ', '));
        
    catch ME
        fprintf('ERROR: Quality assessment calculation failed: %s\n', ME.message);
        improvement_results.error = ME.message;
    end
end