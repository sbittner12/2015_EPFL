function [Preproc] = computePreprocessedEMG(subject, exptype, musclestr, trials, targets)
%COMPUTEMEANMAG Summary of this function goes here
%   Detailed explanation goes here

% Setup coherence calculation parameters

% Maximum index in frequency resolution to show
SUBJECTS = [3,4,5,7,8];
T = .001; % 1 kHz
windowsize = 512;
percentiles = 11;
eventstr = getEventstr(exptype);

[premove_buf, postmove_buf] = getMovementBuffers(SUBJECTS, eventstr);

if isnumeric(subject)
    subject = num2str(subject);
end
ntrials = length(trials);
ntargets = length(targets);
                 
%Setup event string for file I/O
resDir = ['E:',filesep,'Sean', filesep, 'Results', filesep, 'PreprocessedEMG', ... 
          filesep, 'Subject', subject, filesep];
if (exist(resDir, 'dir') ~= 7)
    mkdir(resDir);
end

if (ntrials == 1 && ntargets == 1)
    resfname = [resDir, sprintf('EMGPreproc_%s_%s_trial%d_target%d.mat', exptype, musclestr, trials(1), targets(1))];
elseif (ntrials >= 10 && ntargets == 8)
    resfname = [resDir, sprintf('EMGPreproc_%s_%s_AllTrials.mat', exptype, musclestr)];
elseif ((ntrials == 12) && (ntargets == 1))
    resfname = [resDir, sprintf('EMGPreproc_%s_%s_target%d.mat', exptype, musclestr, targets(1))]; 
else
    fprintf('Error: Compute mean magnitude with all subject trials and targets, or just a single trial.\n');
    resfname = 'dummy.mat';
end

sampling_windows = percentiles;
duration = premove_buf + postmove_buf;
nsamples = round(duration / T)+1;

Preproc_percentile_all = zeros(windowsize, sampling_windows, ntrials*ntargets);
Preproc_uniformlength_all = zeros(nsamples, ntrials*ntargets);

for i=1:ntrials
    trial = trials(i);
    Event = loadEvent(subject, exptype, trial);
    emgdata = extractEMGSignal(subject, exptype, musclestr, trial);
    for j=1:ntargets
        target = targets(j);
        %tic;
        % Extract corresponding eeg and emg signals
        if (ntrials == 1 && ntargets == 1)
            Preproc_reach_duration = extractEMGSignalFast(emgdata, Event, target, 'reach duration');
        end
        
        %emg_uniform_length = extractEMGSignalFast(emgdata, Event, target, ...
%                                                   'uniform length', premove_buf, postmove_buf);
%         Preproc_uniformlength_all(:,(i-1)*ntargets+j) = emg_uniform_length;
        
        emg_percentiles = extractEMGSignalFast(emgdata, Event, target, 'percentiles', windowsize);
        
%         figure, plot(Preproc_reach_duration), title('reach duration');
%         figure, plot(emg_percentiles), title('percentiles');
        
        emgLen = length(emg_percentiles);
        percentile_centers_emg = round(linspace(256, emgLen-256, percentiles));
        percentile_start_inds_emg = percentile_centers_emg - 255;
        percentile_stop_inds_emg = percentile_centers_emg + 256;
        for sampleInd=1:percentiles
            startind_emg = percentile_start_inds_emg(sampleInd);
            stopind_emg = percentile_stop_inds_emg(sampleInd);

            emg_percentile = emg_percentiles(startind_emg:stopind_emg);
            Preproc_percentile_all(:, sampleInd, (i-1)*ntargets+j) = emg_percentile;
        end
    end
end

Preproc_uniform_length = mean(Preproc_uniformlength_all, 2);
Preproc_percentile = mean(Preproc_percentile_all, 3);

if (ntrials == 1 && ntargets == 1)
    save(resfname, 'Preproc_reach_duration', 'Preproc_uniform_length', 'Preproc_percentile'); 
else
    save(resfname, 'Preproc_uniform_length', 'Preproc_percentile'); 
end