function [] = computeSignificanceValues(subject, exptype, musclestr, channels, trials, targets)
%COMPUTECOHERENCE Summary of this function goes here
%   Detailed explanation goes here

% Setup coherence calculation parameters

% Maximum index in frequency resolution to show
T = .001; % 1 kHz
Fs = 1/T;
windowsize = 512;
overlap = 128;
fftbins = (windowsize/2)+1;
pmin = 12;
pmax = 12;
percentiles = 11;
nrandomizations = 19;
blocksize = 100;
% Extraction types can be 'uniform length' or 'reach duration'
extraction_type = 'percentiles'; 
extractionstr = getExtractionstr(extraction_type);

if isnumeric(subject)
    subject = num2str(subject);
end
ntrials = length(trials);
ntargets = length(targets);
                 
%Setup event string for file I/O
resDir = ['E:',filesep,'Sean', filesep, 'Results', filesep, 'Coherence', ... 
          filesep, extractionstr, filesep, 'Subject', subject, filesep, ...
          'SignificanceBootstrapping', filesep];
if (exist(resDir, 'dir') ~= 7)
    mkdir(resDir);
end

if (ntrials == 1 && ntargets == 1)
    resfname = [resDir, sprintf('SignificanceValues_%s_%s_%s_trial%d_target%d.mat', exptype, musclestr, channels, trials(1), targets(1))]; 
elseif (ntrials >= 10 && ntargets == 8)
    resfname = [resDir, sprintf('SignificanceValues_%s_%s_%s_AllTrials.mat', exptype, musclestr, channels)]; 
else
    fprintf('Error: Compute coherence with all subject trials and targets, or just a single trial.\n');
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

CMC_Map = zeros(fftbins, sampling_windows);
% PDC1_Map = zeros(fftbins, sampling_windows);
% PDC2_Map = zeros(fftbins, sampling_windows);
% gPDC1_Map = zeros(fftbins, sampling_windows);
% gPDC2_Map = zeros(fftbins, sampling_windows);
% iCOH1_Map = zeros(fftbins, sampling_windows);
% iCOH2_Map = zeros(fftbins, sampling_windows);

for r=1:nrandomizations
    fprintf('+++++++++++++++++++++ Randomization %d +++++++++++++++++++++\n', r);
    % Calculate coherences for each electrode
    numeegs = getNumEEGs(channels);
    Sxx_all = cell(1,numeegs);
    Syy_all = cell(1,numeegs);
    Sxy_all = cell(1,numeegs);
    % The top row of these cells is 
%     PDC_all = cell(2,numeegs);
%     gPDC_all = cell(2,numeegs);
%     iCOH_all = cell(2,numeegs);
    for j=1:numeegs
        Sxx_all{j} = zeros(fftbins, sampling_windows, ntrials*ntargets);
        Syy_all{j} = zeros(fftbins, sampling_windows, ntrials*ntargets);
        Sxy_all{j} = zeros(fftbins, sampling_windows, ntrials*ntargets);
%         for i=1:2
%             PDC_all{i,j} = zeros(fftbins, sampling_windows, ntrials*ntargets);
%             gPDC_all{i,j} = zeros(fftbins, sampling_windows, ntrials*ntargets);
%             iCOH_all{i,j} = zeros(fftbins, sampling_windows, ntrials*ntargets);
%         end
    end

    EEG = loadEEG(subject, exptype, channels);

    for i=1:ntrials
        trial = trials(i);
        fprintf('Trial %d ', trial);
        tic
        Event = loadEvent(subject, exptype, trial);
        
        eegdata = Phase_randomization(EEG(trial).data');
        %eegdata = EEG(trial).data';
        emgdata = extractEMGSignal(subject, exptype, musclestr, trial);
        srate_new = EEG(trial).srate;
        for j=1:ntargets
            target = targets(j);
            % Extract corresponding eeg and emg signals
            switch extraction_type
                case 'reach duration'
                    eeg = extractEEGSignalFast(eegdata, Event, srate_new, target, extraction_type);
                    emg = extractEMGSignalFast(emgdata, Event, target, extraction_type);
                case 'uniform length'
                    eeg = extractEEGSignalFast(eegdata, Event, srate_new, target, ...
                                           extraction_type, premove_buf, postmove_buf);
                    emg = extractEMGSignalFast(emgdata, Event, target, ...
                                           extraction_type, premove_buf, postmove_buf);
                case 'percentiles'
                    eeg = extractEEGSignalFast(eegdata, Event, srate_new, target, extraction_type, windowsize);
                    emg = extractEMGSignalFast(emgdata, Event, target, extraction_type, windowsize);
%                     nblocks = floor(size(eeg,1)/blocksize);
%                     % trim eeg and emg data
%                     eeg = eeg(1:nblocks*blocksize, :);
%                     emg = emg(1:nblocks*blocksize);
%                     % time block permutation for eeg signal
%                     eeg = Timeblock_permutation(eeg, blocksize);
            end

            numeegs = size(eeg, 2);
            eegLen = length(eeg);
            emgLen = length(emg);
            
            for k=1:numeegs
                if (strcmp(extraction_type, 'reach duration') || strcmp(extraction_type, 'uniform length'))
                    sampleInd = 1;
                    startind = 1;
                    stopind = windowsize;
                    while (stopind < eegLen)
                        x = eeg(startind:stopind, k);
                        y = emg(startind:stopind);

                        % compute directly
                        Sxx = pwelch(x,[],[],windowsize);
                        Syy = pwelch(y,[],[],windowsize);
                        Sxy = cpsd(x,y,[],[],windowsize);

%                         [w, A, C, sbc, fpe, th] = arfit([x, y], pmin, pmax);
%                         acfunc = getacfunc(C,A);
%                         A_w = Autoco2FFT(acfunc, windowsize);

%                         PDC = computePDC(C, A_w, windowsize);
%                         gPDC = computegPDC(C, A_w, windowsize);
%                         iCOH = computeiCOH(C, A_w, windowsize);

                        Sxx_all{k}(:, sampleInd, (i-1)*ntargets+j) = Sxx;
                        Syy_all{k}(:, sampleInd, (i-1)*ntargets+j) = Syy;
                        Sxy_all{k}(:, sampleInd, (i-1)*ntargets+j) = Sxy;
%                         PDC_all{1,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(PDC(1,2,:));
%                         PDC_all{2,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(PDC(2,1,:));
%                         gPDC_all{1,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(gPDC(1,2,:));
%                         gPDC_all{2,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(gPDC(2,1,:));
%                         iCOH_all{1,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(iCOH(1,2,:));
%                         iCOH_all{2,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(iCOH(2,1,:));

                        sampleInd = sampleInd + 1;
                        startind = stopind - overlap + 1;
                        stopind = stopind + windowsize - overlap;
                    end
                elseif (strcmp(extraction_type, 'percentiles'))
                    percentile_centers_eeg = round(linspace(256, eegLen-256, percentiles));
                    percentile_centers_emg = round(linspace(256, emgLen-256, percentiles));
                    percentile_start_inds_eeg = percentile_centers_eeg - 255;
                    percentile_stop_inds_eeg = percentile_centers_eeg + 256;
                    percentile_start_inds_emg = percentile_centers_emg - 255;
                    percentile_stop_inds_emg = percentile_centers_emg + 256;
                    for sampleInd=1:percentiles
                        startind_eeg = percentile_start_inds_eeg(sampleInd);
                        stopind_eeg = percentile_stop_inds_eeg(sampleInd);
                        startind_emg = percentile_start_inds_emg(sampleInd);
                        stopind_emg = percentile_stop_inds_emg(sampleInd);

                        x = eeg(startind_eeg:stopind_eeg, k);
                        y = emg(startind_emg:stopind_emg);
                        
                        % compute directly
                        Sxx = pwelch(x,[],[],windowsize);
                        Syy = pwelch(y,[],[],windowsize);
                        Sxy = cpsd(x,y,[],[],windowsize);

%                         [w, A, C, sbc, fpe, th] = arfit([x, y], pmin, pmax);
%                         acfunc = getacfunc(C,A);
%                         A_w = Autoco2FFT(acfunc, windowsize);
% 
%                         PDC = computePDC(C, A_w, windowsize);
%                         gPDC = computegPDC(C, A_w, windowsize);
%                         iCOH = computeiCOH(C, A_w, windowsize);

                        Sxx_all{k}(:, sampleInd, (i-1)*ntargets+j) = Sxx;
                        Syy_all{k}(:, sampleInd, (i-1)*ntargets+j) = Syy;
                        Sxy_all{k}(:, sampleInd, (i-1)*ntargets+j) = Sxy;
%                         PDC_all{1,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(PDC(1,2,:));
%                         PDC_all{2,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(PDC(2,1,:));
%                         gPDC_all{1,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(gPDC(1,2,:));
%                         gPDC_all{2,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(gPDC(2,1,:));
%                         iCOH_all{1,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(iCOH(1,2,:));
%                         iCOH_all{2,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(iCOH(2,1,:));
                    end
                end
            end
        end
        fprintf('took %f seconds\n', toc);
    end

    % Take the cross trial average and compute corticomuscular coherence
    for k=1:numeegs
        Sxx_mean_k = mean(Sxx_all{k}, 3);
        Syy_mean_k = mean(Syy_all{k}, 3);
        Sxy_mean_k = mean(Sxy_all{k}, 3);
        CMC_mean{k} = (abs(Sxy_mean_k).^2) ./ (Sxx_mean_k.*Syy_mean_k);

%         PDC_mean{1,k} = mean(PDC_all{1,k}, 3);
%         PDC_mean{2,k} = mean(PDC_all{2,k}, 3);
%         gPDC_mean{1,k} = mean(gPDC_all{1,k}, 3);
%         gPDC_mean{2,k} = mean(gPDC_all{2,k}, 3);
%         iCOH_mean{1,k} = mean(iCOH_all{1,k}, 3);
%         iCOH_mean{2,k} = mean(iCOH_all{2,k}, 3);
    end

    % Generate dynamic coherence maps by comparing against a significance
    % threshold th
    temporal_samples = size(CMC_mean{1}, 2);
    if (temporal_samples ~= sampling_windows)
        fprintf('Incorrect sampling window calcuation.\n')
        return;
    end
    CMC_mat = zeros(fftbins, sampling_windows, numeegs);

%     PDC1_mat = zeros(fftbins, sampling_windows, numeegs);
%     PDC2_mat = zeros(fftbins, sampling_windows, numeegs);
%     gPDC1_mat = zeros(fftbins, sampling_windows, numeegs);
%     gPDC2_mat = zeros(fftbins, sampling_windows, numeegs);
%     iCOH1_mat = zeros(fftbins, sampling_windows, numeegs);
%     iCOH2_mat = zeros(fftbins, sampling_windows, numeegs);

    for k=1:numeegs
        CMC_mat(:,:,k) = CMC_mean{k};
%         PDC1_mat(:,:,k) = PDC_mean{1,k};
%         PDC2_mat(:,:,k) = PDC_mean{2,k};
%         gPDC1_mat(:,:,k) = gPDC_mean{1,k};
%         gPDC2_mat(:,:,k) = gPDC_mean{2,k};
%         iCOH1_mat(:,:,k) = iCOH_mean{1,k};
%         iCOH2_mat(:,:,k) = iCOH_mean{2,k};
    end

    % Average the results across electrodes
    CMC_Map = CMC_Map + mean(CMC_mat, 3);
%     PDC1_Map = PDC1_Map + mean(PDC1_mat, 3);
%     PDC2_Map = PDC2_Map + mean(PDC2_mat, 3);
%     gPDC1_Map = gPDC1_Map + mean(gPDC1_mat, 3);
%     gPDC2_Map = gPDC2_Map + mean(gPDC2_mat, 3);
%     iCOH1_Map = iCOH1_Map + mean(iCOH1_mat, 3);
%     iCOH2_Map = iCOH2_Map + mean(iCOH2_mat, 3);
end

CMC_Sig = CMC_Map / nrandomizations;
% PDC1_Sig = PDC1_Map / nrandomizations;
% x1 = PDC1_Sig(1,1);
% x2 = PDC1_Sig(2,1);
% x3 = PDC1_Sig(3,1);
% x4 = PDC1_Sig(4,1);
% fprintf('Final %d: %f %f %f %f \n', r, x1, x2, x3, x4);
% PDC2_Sig = PDC2_Map / nrandomizations;
% gPDC1_Sig = gPDC1_Map / nrandomizations;
% gPDC2_Sig = gPDC2_Map / nrandomizations;
% iCOH1_Sig = iCOH1_Map / nrandomizations;
% iCOH2_Sig = iCOH2_Map / nrandomizations;

save(resfname, 'CMC_Sig'); %, 'PDC1_Sig', 'PDC2_Sig', 'gPDC1_Sig', 'gPDC2_Sig', 'iCOH1_Sig', 'iCOH2_Sig');