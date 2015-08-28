function [Mean_Mag] = computeMeanMag(subject, exptype, musclestr, trials, targets)
%COMPUTEMEANMAG Summary of this function goes here
%   Detailed explanation goes here

% Setup coherence calculation parameters

% Maximum index in frequency resolution to show
T = .001; % 1 kHz
windowsize = 512;
percentiles = 11;
% Extraction types can be 'uniform length' or 'reach duration'
extraction_type = 'percentiles'; 
extractionstr = getExtractionstr(extraction_type);

if isnumeric(subject)
    subject = num2str(subject);
end
ntrials = length(trials);
ntargets = length(targets);
                 
%Setup event string for file I/O
resDir = ['E:',filesep,'Sean', filesep, 'Results', filesep, 'MeanMag', ... 
          filesep, extractionstr, filesep, 'Subject', subject, filesep];
if (exist(resDir, 'dir') ~= 7)
    mkdir(resDir);
end

if (ntrials == 1 && ntargets == 1)
    resfname = [resDir, sprintf('MeanMag_%s_%s_%s_trial%d_target%d.mat', exptype, musclestr, trials(1), targets(1))]; 
elseif (ntrials >= 10 && ntargets == 8)
    resfname = [resDir, sprintf('MeanMag_%s_%s_AllTrials.mat', exptype, musclestr)]; 
else
    fprintf('Error: Compute mean magnitude with all subject trials and targets, or just a single trial.\n');
    resfname = 'dummy.mat';
end

eventstr = getEventstr(exptype);

% Initialize sampling window properties
switch extraction_type
    case 'reach duration'
        % Set to max samples
        % End of signals past 14 don't have proper meaning
        % Use uniform length to compare
        sampling_windows = 22; 
    case 'uniform length'
        [premove_buf, postmove_buf] = getMovementBuffers(SUBJECTS, eventstr);
        fprintf('Using uniform length for computing corticomuscular coherence.\n')
        fprintf('Maximum pre movement onset buffer: %.2fs\n', premove_buf);
        fprintf('Maximum post movement onset buffer: %.2fs\n', postmove_buf);
        duration = premove_buf + postmove_buf;
        nsamples = round(duration / T);
        assert(nsamples > 511);
        sampling_windows = floor((nsamples-512) / (windowsize-overlap)) + 1;
        samples_left_out =  mod((nsamples-512), (windowsize-overlap));
        fprintf('Sampling windows: %d\n', sampling_windows);
        fprintf('Samples left out: %d\n', samples_left_out);
        postmove_buf_actual = postmove_buf - samples_left_out*T;
        % Must change due to imperfect sampling window size
        fprintf('Actual post movement onset buffer: %.2fs\n', postmove_buf_actual);
    case 'percentiles'
        sampling_windows = percentiles;
    otherwise
        fprintf('Error: extraction type "%s" invalid.\n', extraction_type);
        return;
        
end

Mag_all = zeros(windowsize, sampling_windows, ntrials*ntargets);
MeanMag_all = zeros(windowsize, sampling_windows, ntrials*ntargets);

for i=1:ntrials
    trial = trials(i);
    fprintf('Trial %d   ', trial);
    tic;
    Event = loadEvent(subject, exptype, trial);
    emgdata = extractEMGSignal(subject, exptype, musclestr, trial);
    for j=1:ntargets
        %tic;
        target = targets(j);
        % Extract corresponding eeg and emg signals
        switch extraction_type
            case 'reach duration'
                emg = extractEMGSignal(emgdata, Event, target, extraction_type);
            case 'uniform length'
                emg = extractEMGSignalFast(emgdata, Event, target, ...
                                       extraction_type, premove_buf, postmove_buf);
            case 'percentiles'
                emg = extractEMGSignalFast(emgdata, Event, target, extraction_type, windowsize);
        end
        
        emgLen = length(emg);
        if (strcmp(extraction_type, 'percentiles'))
            percentile_centers_emg = round(linspace(256, emgLen-256, percentiles));
            percentile_start_inds_emg = percentile_centers_emg - 255;
            percentile_stop_inds_emg = percentile_centers_emg + 256;
            for sampleInd=1:percentiles
                startind_emg = percentile_start_inds_emg(sampleInd);
                stopind_emg = percentile_stop_inds_emg(sampleInd);
                
                emg_percentile = emg(startind_emg:stopind_emg);
                Mag_all(:, sampleInd, (i-1)*ntargets+j) = emg_percentile;
                MeanMag_all(:,sampleInd,(i-1)*ntargets+j) = est_amp(emg_percentile, 2, 60); % paremeters used in mbPLS
            end
        end
    fprintf(' %.2fs\n', toc);
    end
end

Mean_Mag = mean(MeanMag_all, 3);

switch extraction_type
    case 'reach duration'
        save(resfname, 'CMC_Map', 'PDC1_Map', 'PDC2_Map', 'gPDC1_Map', 'gPDC2_Map', 'iCOH1_Map', 'iCOH2_Map');
    case 'uniform length'
        save(resfname, 'CMC_Map', 'PDC1_Map', 'PDC2_Map', 'gPDC1_Map', 'gPDC2_Map', 'iCOH1_Map', 'iCOH2_Map', 'premove_buf', 'postmove_buf_actual');
    case 'percentiles'
        save(resfname, 'Mean_Mag'); %, 'PDC1_Map', 'PDC2_Map', 'gPDC1_Map', 'gPDC2_Map', 'iCOH1_Map', 'iCOH2_Map');
end
