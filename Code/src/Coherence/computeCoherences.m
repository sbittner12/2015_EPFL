function [CMC] = computeCoherences(subject, exptype, musclestr, channels, trials, targets)
%COMPUTECOHERENCE Summary of this function goes here
%   Detailed explanation goes here

% Setup coherence calculation parameters

% Subject that will be used to determine pre movement onset buffers and 
% post movement onset buffers
SUBJECTS = [3,4,7];
% Maximum index in frequency resolution to show
MAX_IND_SHOW = 26;  % shows up to 50.96 Hz
FREQ_TO_PLOT = [10, 20, 40, 50];
T = .001; % 1 kHz
srate = 1000;
windowsize = 512;
overlap = 128;
fftbins = (windowsize/2)+1;
freq_res = srate/windowsize;
EEG_T = 1 / 2048; 
pmin = 12;
pmax = 12;
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
resDir = ['E:',filesep,'Sean', filesep, 'Results', filesep, 'Coherence', ... 
          filesep, extractionstr, filesep, 'Subject', subject, filesep, ...
          'RawCoherence', filesep, musclestr, filesep];
if (exist(resDir, 'dir') ~= 7)
    mkdir(resDir);
end

if (ntrials == 1 && ntargets == 1)
    resfname = [resDir, sprintf('Coherences_%s_%s_trial%d_target%d.mat', exptype, channels, trials(1), targets(1))]; 
elseif (ntrials >= 10 && ntargets == 8)
    resfname = [resDir, sprintf('Coherences_%s_%s_AllTrials1.mat', exptype, channels)]; 
elseif ((ntrials == 12) && (ntargets == 1))
    resfname = [resDir, sprintf('Coherences_%s_%s_target%d.mat', exptype, channels, targets(1))]; 
else
    fprintf('Error: Compute CMC with all subject trials and targets, or just a single trial.\n');
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

% Calculate coherences for each electrode
numeegs = getNumEEGs(channels);
Sxx_all = cell(1,numeegs);
Syy_all = cell(1,numeegs);
Sxy_all = cell(1,numeegs);
% The top row of these cells is 
% PDC_all = cell(2,numeegs);
% gPDC_all = cell(2,numeegs);
% iCOH_all = cell(2,numeegs);
for j=1:numeegs
    Sxx_all{j} = zeros(fftbins, sampling_windows, ntrials*ntargets);
    Syy_all{j} = zeros(fftbins, sampling_windows, ntrials*ntargets);
    Sxy_all{j} = zeros(fftbins, sampling_windows, ntrials*ntargets);
%     for i=1:2
%         PDC_all{i,j} = zeros(fftbins, sampling_windows, ntrials*ntargets);
%         gPDC_all{i,j} = zeros(fftbins, sampling_windows, ntrials*ntargets);
%         iCOH_all{i,j} = zeros(fftbins, sampling_windows, ntrials*ntargets);
%     end
end


EEG = loadEEG(subject, exptype, channels, srate);

for i=1:ntrials
    trial = trials(i);
%     fprintf('Trial %d   ', trial);
%     tic;
    Event = loadEvent(subject, exptype, trial);
    eegdata = EEG(trial).data';
    emgdata = extractEMGSignal(subject, exptype, musclestr, trial);
    srate_new = EEG(trial).srate;
        for j=1:ntargets
            %tic;
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
%                     fprintf('stopind: %d\n', stopind);
%                     fprintf('eeglen: %d\n', eegLen);
%                     fprintf('emgLen: %d\n', emgLen);
                    x = eeg(startind:stopind, k);
                    y = emg(startind:stopind);

                    % compute directly
                    Sxx = pwelch(x,[],[],windowsize);
                    Syy = pwelch(y,[],[],windowsize);
                    Sxy = cpsd(x,y,[],[],windowsize);
% 
%                     [w, A, C, sbc, fpe, th] = arfit([x', y'], pmin, pmax);
%                     acfunc = getacfunc(C,A);
%                     A_w = Autoco2FFT(acfunc, windowsize);
% 
%                     PDC = computePDC(C, A_w, windowsize);
%                     gPDC = computegPDC(C, A_w, windowsize);
%                     iCOH = computeiCOH(C, A_w, windowsize);

                    Sxx_all{k}(:, sampleInd, (i-1)*ntargets+j) = Sxx;
                    Syy_all{k}(:, sampleInd, (i-1)*ntargets+j) = Syy;
                    Sxy_all{k}(:, sampleInd, (i-1)*ntargets+j) = Sxy;
%                     PDC_all{1,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(PDC(1,2,:));
%                     PDC_all{2,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(PDC(2,1,:));
%                     gPDC_all{1,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(gPDC(1,2,:));
%                     gPDC_all{2,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(gPDC(2,1,:));
%                     iCOH_all{1,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(iCOH(1,2,:));
%                     iCOH_all{2,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(iCOH(2,1,:));

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
                    
                    if (sampleInd == 14)
                        figure, plot(x);
                        figure, plot(y);
                        figure, plot(Sxx)
                        figure, plot(Syy)
                        figure, plot(1:257, abs(Sxy))
                        CMC_mean = (abs(Sxy).^2) ./ (abs(Sxx).*abs(Syy));
                        figure, imagesc(CMC_mean)
                        title('sample 6');
                    end
%                     [w, A, C, sbc, fpe, th] = arfit([x, y], pmin, pmax);
%                     acfunc = getacfunc(C,A);
%                     A_w = Autoco2FFT(acfunc, windowsize);
% 
%                     PDC = computePDC(C, A_w, windowsize);
%                     gPDC = computegPDC(C, A_w, windowsize);
%                     iCOH = computeiCOH(C, A_w, windowsize);

                    Sxx_all{k}(:, sampleInd, (i-1)*ntargets+j) = Sxx;
                    Syy_all{k}(:, sampleInd, (i-1)*ntargets+j) = Syy;
                    Sxy_all{k}(:, sampleInd, (i-1)*ntargets+j) = Sxy;
%                     PDC_all{1,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(PDC(1,2,:));
%                     PDC_all{2,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(PDC(2,1,:));
%                     gPDC_all{1,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(gPDC(1,2,:));
%                     gPDC_all{2,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(gPDC(2,1,:));
%                     iCOH_all{1,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(iCOH(1,2,:));
%                     iCOH_all{2,k}(:, sampleInd, (i-1)*ntargets+j) = squeeze(iCOH(2,1,:));
                end
            end
%             CMC_mean = (abs(Sxy_all{k}).^2) ./ (abs(Sxx_all{k}).*abs(Syy_all{k}));
%             figure, imagesc(CMC_mean)
%             title('electrode mean');
        end
        end
%     fprintf(' %.2fs\n', toc);
end

% Take the cross trial average and compute corticomuscular coherence
for k=1:numeegs
    Sxx_mean_k = mean(Sxx_all{k}, 3);
    Syy_mean_k = mean(Syy_all{k}, 3);
    Sxy_mean_k = mean(Sxy_all{k}, 3);
    CMC_mean{k} = (abs(Sxy_mean_k).^2) ./ (Sxx_mean_k.*Syy_mean_k);
    
%     PDC_mean{1,k} = mean(PDC_all{1,k}, 3);
%     x1 = PDC_mean{1,k}(1,1);
%     x2 = PDC_mean{1,k}(2,1);
%     x3 = PDC_mean{1,k}(3,1);
%     x4 = PDC_mean{1,k}(4,1);
%     fprintf('electrode %d: %f %f %f %f   mean\n', k, x1, x2, x3, x4);
%     PDC_mean{2,k} = mean(PDC_all{2,k}, 3);
%     gPDC_mean{1,k} = mean(gPDC_all{1,k}, 3);
%     gPDC_mean{2,k} = mean(gPDC_all{2,k}, 3);
%     iCOH_mean{1,k} = mean(iCOH_all{1,k}, 3);
%     iCOH_mean{2,k} = mean(iCOH_all{2,k}, 3);
end

% Generate dynamic coherence maps by comparing against a significance
% threshold th
temporal_samples = size(CMC_mean{1}, 2);
if (temporal_samples ~= sampling_windows)
    fprintf('Incorrect sampling window calcuation.\n')
    return;
end
CMC_mat = zeros(fftbins, sampling_windows, numeegs);
% size(CMC_mat)
% size(CMC_mean)
% size(CMC_mean{1})
% PDC1_mat = zeros(fftbins, sampling_windows, numeegs);
% PDC2_mat = zeros(fftbins, sampling_windows, numeegs);
% gPDC1_mat = zeros(fftbins, sampling_windows, numeegs);
% gPDC2_mat = zeros(fftbins, sampling_windows, numeegs);
% iCOH1_mat = zeros(fftbins, sampling_windows, numeegs);
% iCOH2_mat = zeros(fftbins, sampling_windows, numeegs);

for k=1:numeegs
    CMC_mat(:,:,k) = CMC_mean{k};
%     PDC1_mat(:,:,k) = PDC_mean{1,k};
%     PDC2_mat(:,:,k) = PDC_mean{2,k};
%     gPDC1_mat(:,:,k) = gPDC_mean{1,k};
%     gPDC2_mat(:,:,k) = gPDC_mean{2,k};
%     iCOH1_mat(:,:,k) = iCOH_mean{1,k};
%     iCOH2_mat(:,:,k) = iCOH_mean{2,k};
end

% Average the results across electrodes
CMC_Map = mean(CMC_mat, 3);
% PDC1_Map = mean(PDC1_mat, 3);
% x1 = PDC1_Map(1,1);
% x2 = PDC1_Map(2,1);
% x3 = PDC1_Map(3,1);
% x4 = PDC1_Map(4,1);
% fprintf('average: %f %f %f %f \n', x1, x2, x3, x4);
% PDC2_Map = mean(PDC2_mat, 3);
% gPDC1_Map = mean(gPDC1_mat, 3);
% gPDC2_Map = mean(gPDC2_mat, 3);
% iCOH1_Map = mean(iCOH1_mat, 3);
% iCOH2_Map = mean(iCOH2_mat, 3);

switch extraction_type
    case 'reach duration'
        save(resfname, 'CMC_Map'); %, 'PDC1_Map', 'PDC2_Map', 'gPDC1_Map', 'gPDC2_Map', 'iCOH1_Map', 'iCOH2_Map');
    case 'uniform length'
        save(resfname, 'CMC_Map', 'PDC1_Map', 'PDC2_Map', 'gPDC1_Map', 'gPDC2_Map', 'iCOH1_Map', 'iCOH2_Map', 'premove_buf', 'postmove_buf_actual');
    case 'percentiles'
        save(resfname, 'CMC_Map'); %, 'PDC1_Map', 'PDC2_Map', 'gPDC1_Map', 'gPDC2_Map', 'iCOH1_Map', 'iCOH2_Map');
end

% Maps = {CMC_Map, PDC1_Map, PDC2_Map, gPDC1_Map, gPDC2_Map, iCOH1_Map, iCOH2_Map};
% mapnames = {'Corticomuscular coherence (CMC)', ...
%             'Generalized PDC (GPDC) EEG -> EMG', 'Generalized PDC (GPDC) EMG -> EEG', ...
%             'Partial directed coherence EEG -> EMG (PDC)', 'Partial directed coherence EMG -> EEG (PDC)', ...
%             'independent coherence (iCOH) EEG -> EMG', 'independent coherence (iCOH) EMG -> EEG'};
% N = length(Maps);
% for i=1:N
%     map = Maps{i};
%     mapname = mapnames{i};
% 
%     figure, imagesc(map(1:MAX_IND_SHOW,:));
%     title(mapname);
%     set(gca,'YDir','normal');
%     %% Setup axes
%     % Y axis
%     ax = gca;
%     yaxis_pos = (FREQ_TO_PLOT / freq_res);
%     set(ax,'YTick',yaxis_pos);
%     % make tick labels
%     nlabels = length(FREQ_TO_PLOT);
%     yaxis_labels = cell(1,nlabels); 
%     for i=1:nlabels
%         yaxis_labels{i} = sprintf('%d Hz', FREQ_TO_PLOT(i));
%     end
%     set(ax,'YTickLabel',yaxis_labels);
%     hold on;
% 
%     % Plot the start and uniform length buffer sizes
%     if (strcmp(extraction_type, 'uniform length'))
%         duration = premove_buf + postmove_buf_actual;
%         bordershiftl = .5; % used to place trigger at beginning of plot
%         bordershiftr = .5; % used to place end2 at end
%         begin_pos = 1-bordershiftl;
%         end_pos = sampling_windows+bordershiftr;
%         start_pos = sampling_windows*(premove_buf/duration);
%         % X axis
%         set(ax, 'XTick', [begin_pos, start_pos, end_pos]);
%         xaxis_labels{1} = sprintf('-%.2fs', premove_buf);
%         xaxis_labels{2} = sprintf('move onset');
%         xaxis_labels{3} = sprintf('+%.2fs', postmove_buf_actual);
%         set(ax, 'XTickLabel',xaxis_labels);
% 
%         linewidth = 3;
%         green_color = [34,139,34]/255;
%         plot([start_pos, start_pos], [0, MAX_IND_SHOW+1], 'Color', green_color, 'LineWidth', linewidth);
%     % Plot the trial trigger, start, and end points    
%     elseif (strcmp(extraction_type, 'reach duration') && (ntrials == 1) && (ntargets == 1))
%         % X axis
%         set(ax, 'XTick', []);
%         dataDir = ['E:',filesep,'Sean', filesep, 'Data', filesep];
%         subjectEventDir = [dataDir, 'Event', filesep, 'Subject_', subject, filesep]; 
%         if (strcmp(exptype, 'Healthy'))
%             eventstr = 'HM';
%         else
%             eventstr = exptype;
%         end
%         eventfname = [subjectEventDir, sprintf('%s_%d_event', eventstr, trial)];
%         load(eventfname);
% 
%         % Get the right triggers since sometimes there are too many.
%         trigger = parseTriggers(Event.EEG.Trigger, Event.EEG.End(end));
%         trigger_ind = trigger(target);
%         start1_ind = Event.EEG.Start((2*target)-1);
%         end1_ind   = Event.EEG.End((2*target)-1);
%         start2_ind = Event.EEG.Start(2*target);
%         end2_ind =   Event.EEG.End(2*target);
% 
%         trigger_t  = EEG_T * trigger_ind;
%         start1_t   = EEG_T * start1_ind;
%         end1_t   = EEG_T * end1_ind;
%         start2_t   = EEG_T * start2_ind;
%         end2_t   = EEG_T * end2_ind;
%         duration = end2_t - trigger_t;
% 
%         bordershiftl = .4; % used to place trigger at beginning of plot
%         bordershiftr = .45; % used to place end2 at end
%         trig_pos = 1-bordershiftl;
%         start1_pos = temporalSamples*((start1_t - trigger_t) / duration);
%         end1_pos = temporalSamples*((end1_t - trigger_t) / duration);
%         start2_pos = temporalSamples*((start2_t - trigger_t) / duration);
%         end2_pos = temporalSamples+bordershiftr;
%         overlapshift = .05;
%         if (abs(end1_pos - start2_pos) < .1)
%             end1_pos = end1_pos - overlapshift;
%             start2_pos = start2_pos + overlapshift;
%         end
%         if (MAX_IND_SHOW == 26)
%             textshiftv = .5;
%         elseif (MAX_IND_SHOW == 257)
%             textshiftv = 5;
%         end
% 
%         linewidth = 3;
%         textshifth = .5;
%         green_color = [34,139,34]/255;
%         plot([trig_pos, trig_pos], [0, MAX_IND_SHOW+1], 'b', 'LineWidth', linewidth);
%         text(trig_pos - textshifth, -textshiftv, sprintf('Trigger: %.2fs', trigger_t), 'Color', 'b');
% 
%         plot([start1_pos, start1_pos], [0, MAX_IND_SHOW+1], 'Color', green_color, 'LineWidth', linewidth);
%         text(start1_pos - textshifth, -textshiftv, sprintf('Start: %.2fs', start1_t), 'Color', green_color);
% 
%         plot([end1_pos, end1_pos], [0, MAX_IND_SHOW+1], 'r', 'LineWidth', linewidth);
%         text(end1_pos - textshifth, MAX_IND_SHOW+2*textshiftv, sprintf('Stop: %.2fs', end1_t), 'Color', 'r');
% 
%         plot([start2_pos, start2_pos], [0, MAX_IND_SHOW+1], 'Color', green_color, 'LineWidth', linewidth);
%         text(start2_pos - textshifth, -textshiftv, sprintf('Start: %.2fs', start2_t), 'Color', green_color);
% 
%         plot([end2_pos, end2_pos], [0, MAX_IND_SHOW+1], 'r', 'LineWidth', linewidth);
%         text(end2_pos - textshifth, MAX_IND_SHOW+2*textshiftv, sprintf('Stop: %.2fs', end2_t), 'Color', 'r');
%     end
% 
% end
% 
