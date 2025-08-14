# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

WaS4 is a MATLAB-based data processing pipeline for multimodal physiological research combining EEG, eye-tracking (Pupil Labs), motion capture (XSens), and cardiac data (ECG). The repository contains three sequential processing scripts that align and merge heterogeneous sensor data streams.

## Architecture

### Data Processing Pipeline
The system follows a streamlined two-stage processing architecture:

1. **WaS4_01_xSensToMat.m**: Conversion of XSens Excel data to MATLAB format (preprocessing step)
2. **WaS4_00_pupilMerge.m**: Primary alignment and merging of EEG, Pupil Labs, and XSens data (integrated)

### Data Flow
```
Raw Data Sources → XSens Conversion → Integrated Processing → Final Dataset
- EEG (.bdf)         ↘                     ↗
- XDF streams         ↘                   ↗
- Pupil Labs (.csv)    → pupilMerge → Complete EEG structure  
- XSens (.xlsx) → xSensToMat ↗              (with all modalities)
```

**XSens Integration Details:**
- **Cross-correlation alignment**: Uses `xcorr()` to align .mat file data with XSens stream timestamps  
- **Real timestamps**: Leverages actual XSens stream timestamps from XDF for precise temporal alignment
- **Signal correlation**: Aligns Y-axis CoM data between stream and .mat file for optimal synchronization
- **Size validation**: Automatic padding/truncation for sample count matching with EEG data
- **Robust fallback**: Linear interpolation backup when stream data unavailable
- **Channels added**: Center of Mass position/velocity, foot position, and toe position data

### Time Synchronization Strategy
The pipeline uses a sophisticated multi-method time alignment approach:
- **Primary**: LSL (Lab Streaming Layer) sync events with sequence numbers
- **Fallback**: recording.begin/end events
- **Quality Control**: Cross-correlation validation, RMSE < 0.1s, R² > 0.99

### File Organization
- Subject folders follow pattern `WaS_XXX` where XXX is zero-padded subject number
- Raw data in `/DATA/` directory structure
- Aligned output in dedicated output directories
- Processing maintains full audit trail with quality metrics

## Development Commands

### Running the Pipeline
Execute scripts in sequence:
```matlab
% Stage 1: Convert XSens data to MATLAB format
WaS4_01_xSensToMat

% Stage 2: Integrated alignment and merging pipeline
WaS4_00_pupilMerge  % Now includes XSens integration
```

### Optional Additional Analysis
```matlab
% Gait event detection (if XSens foot data available)
WaS4_03_gaitAnalysis
```

### Key Configuration Variables
Each script contains configuration sections at the top:
- `dataPath`: Input data directory
- `outputPath`/`alignedPath`/`outPath`: Output directories
- `processAllSubjects`: Boolean for batch vs. selective processing
- `specificSubjects`: Array for selective processing
- `overwriteExisting`: Boolean for reprocessing
- `verboseOutput`: Boolean for detailed logging

### Quality Control Parameters
- `minSyncEvents`: Minimum sync events required (typically 2)
- `maxRMSE`: Maximum time alignment error (0.1-0.5s)
- `minRSquared`: Minimum correlation for good alignment (0.95-0.99)

## Key Functions and Utilities

### Time Synchronization Functions
- `match_sync_events_robust()`: Primary LSL sync event matching
- `match_recording_events()`: Fallback recording event matching
- `csv_to_xdf_time()`: Time domain conversion functions

### Data Loading Functions
- `read_large_excel_safe()`: Enhanced Excel reader with network drive optimization
- `find_timestamp_column()`: Flexible timestamp column detection
- `extract_special_events()`: Experimental event filtering

### Quality Assessment
- `create_alignment_plots()`: Comprehensive alignment quality visualization
- `create_merge_plots()`: Data integration quality assessment

## Error Handling Strategy

The pipeline implements robust error handling:
- Graceful degradation when data sources are missing
- Comprehensive logging with detailed error reporting
- Continuation processing (failed subjects don't stop batch processing)
- Automatic fallback methods for critical operations

## Data Quality Validation

### Temporal Alignment Validation
- Multi-method sync event matching with sequence number validation
- Statistical quality metrics (RMSE, R-squared) for each subject
- Cross-correlation validation for ECG alignment
- Interpolation coverage analysis for merged data channels

### File Handling Optimizations
- Local copy strategy for large files on network drives
- Multiple Excel reading methods with automatic fallback
- Robust file system error handling with cleanup procedures

## Dependencies

### Required MATLAB Toolboxes
- EEGLAB (EEG data processing)
- Signal Processing Toolbox (interpolation, filtering)
- Statistics and Machine Learning Toolbox (correlation analysis)

### External Libraries
- `load_xdf`: XDF file reading capability
- `pop_biosig`: EEG file format support (Biosig)
- `edfread`/`edfinfo`: EDF file reading for ECG data

### Data Format Requirements
- EEG: BioSemi BDF format with counter channel ('CNT')
- LSL: XDF format with marker/event streams
- Pupil Labs: CSV format with nanosecond timestamps
- XSens: Excel format with specific sheet structure
- ECG: EDF format from Faros devices

## Performance Considerations

### Large File Handling
- Scripts implement local copy strategy for network drives
- Progress reporting for large data processing
- Memory-efficient processing with streaming for large datasets
- Compressed output format (-v7.3) for large merged datasets

### Processing Optimization
- Parallel tool calls where possible
- Efficient interpolation methods (cubic, makima)
- Selective channel processing to minimize memory usage