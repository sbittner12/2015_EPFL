% function [FS_trft,FO_trft] = Analysis_rat_608_SNR()

clear all;

subjects = 3:9;

for subject = subjects
    fprintf('Calculating SNRs for subject %d\n', subject);
    % Directory that holds the data from the session
    useBlocks = 1:12;

    plotDiagnostics = true;
    exptype = 'OAF';
    H.noOfBlocks = length(useBlocks);
    H.datasetName = 'Subject 3 outside exoskeleton';
    H.defectiveEl = 0;
    H.noOfEl = 46;
    H.CARch = [];
    H.electrodeConfig = {1, 2, 3, 4, 5, 6; ...
                         7, 8, 9, 10, 11, 12; ...
                         13, 14, 15, 16, 17, 18; ...
                         19, 20, 21, 22, 23, 24; ...
                         25, 26, 27, 28, 29, 30; ...
                         31, 32, 33, 34, 35, 36; ...
                         37, 38, 39, 40, 41, 42; ...
                         43, 44, 45, 46, [], []};

    electrodeList = 1:H.noOfEl;
    H.sampleRate = 1000;
    H.sampleRate_old = 2048;

    %% Collecting raw data file names
    % These files contain 30kHz broadband recordings (ns5)
    EEG.data = cell(length(useBlocks),1);
    trig.move = cell(length(useBlocks),1);
    trig.calm = cell(length(useBlocks),1);
    trig.task = cell(length(useBlocks),1);
    X = loadEEG(subject, exptype, 'all');
    for blk = 1:length(useBlocks)
        EEG.data{blk} = X(blk).data';
        Event = loadEvent(subject, exptype, blk);
        nepochs = length(Event.EEG.Start) / 2;
        triggers = round(parseTriggers(Event.EEG.Trigger, Event.EEG.End(16)) * H.sampleRate / H.sampleRate_old);
        startinds = round(Event.EEG.Start  * H.sampleRate / H.sampleRate_old);
        endinds = round(Event.EEG.End  * H.sampleRate / H.sampleRate_old);
        % Collect indices for baseline
        baselineInds{blk} = [1, triggers(1);
                            endinds(2:2:14)', triggers(2:8)';
                            endinds(16), size(EEG.data{blk}, 1)];
        % Collect indices for forward reaching epochs
        trig.SF{blk} = [triggers.', endinds(1:2:15)'];
        % Collect indices for backward reaching epochs
        trig.SB{blk} = [startinds(2:2:16)', endinds(2:2:16)'];
    end



    %% Doing the STFT analysis

    trialPast = round(0.1 * H.sampleRate);
    trialFuture = round(1.0 * H.sampleRate);
    trialTime = (-trialPast:trialFuture)/H.sampleRate;
    fftWinSize = 512; %ceil(H.sampleRate/10);
    fftWinFunction = hamming(fftWinSize);
    fftStep = fftWinSize / 4; %floor(H.sampleRate/40);
    fftBaseStep = floor(fftWinSize/2);
    fftTrialPast = trialPast + floor(fftWinSize/2);
    fftTrialFuture = trialFuture + floor(fftWinSize/2);
    offsetBefore = fftTrialPast;
    offsetAfter = fftTrialFuture;

    usedElectrodes = 1:H.noOfEl;

    %%

    baseOffset = ceil(H.sampleRate/10);

    [trft_base,trft_frequencies,~,trft_blkInd] = baselineAmplitude_SE(EEG.data, ...
                                                                        [], ...
                                                                        baselineInds, ...   
                                                                        baseOffset, ...
                                                                        fftWinSize, ...
                                                                        fftBaseStep, ...
                                                                        H.sampleRate, ...
                                                                        fftWinFunction);

    % trft_base is unnormalized  (Nf x Nc x Nt)
    trft_baseMean = mean(trft_base,3);   % normalized across trials
    trft_baseStd = zeros(size(trft_baseMean));
    for ch = 1:length(electrodeList)
        trft_baseStd(:,ch) = std(trft_base(:,ch,:),0,3);   % standard dev across trials
    end

    trft_baseNormStd = trft_baseStd./trft_baseMean;

    %%

    freqLim = 50;
    freqLimInd = find(trft_frequencies < freqLim,1,'last');

    clear param;
    param.xAxis = trft_frequencies(1:freqLimInd);
    param.design = H.electrodeConfig;
    param.normalize = 'false';
    param.pictureTitle = 'Normalized standard deviation for t7.2014.03.19';
    param.xlabel = 'Frequency (Hz)';
    param.ylabel = 'Amplitude standard deviation';
    param.markedCh = H.CARch;
    %plotTimeMultiarray({trft_baseNormStd(1:freqLimInd,:)},param);

    %% Plotting the spectrogram

    freqLim = 50;
    freqLimInd = find(trft_frequencies < freqLim,1,'last');

    clear param;
    param.normalize = 'false';
    param.xAxis = trft_frequencies(1:freqLimInd);
    param.design = H.electrodeConfig;
    param.markedCh = H.CARch;
    %plotTimeMultiarray({log(trft_baseMean(1:freqLimInd,:))},param);

    %% Plotting the trft of the whole block together with the cues


    freqLim = 51;
    freqLimInd = find(trft_frequencies < freqLim,1,'last');

    pltStep = floor(fftWinSize/16);
    % meanEl = setdiff(1:H.noOfEl,H.noisyEl);

    eeggroups = {'frontal', 'parietal', 'left', 'right', 'center'};
    Neegs = length(eeggroups);
    for i=1:Neegs
        eeggroup = eeggroups{i};
        chaninds = getEEGChannelInds(X(1).chanlocs, eeggroup);
        meanEl = 1:length(chaninds);
        for j=1:H.noOfBlocks
            eegdata{j} = EEG.data{j}(:,chaninds);
        end
        for blk = 1:1
            [trftOfBlock,~,trftOfBlock_time] = baselineAmplitude_SE(eegdata, ...
                                                                trft_baseMean, ... 
                                                                [], ... 
                                                                baseOffset, ...
                                                                fftWinSize, ...
                                                                pltStep, ...
                                                                H.sampleRate, ...
                                                                fftWinFunction, ...
                                                                blk);

            trftOfBlockMean = squeeze(mean(trftOfBlock(:,meanEl,:),2));
            continue
            figure
            hold on
            imagesc(trftOfBlock_time/H.sampleRate, trft_frequencies,trftOfBlockMean)
            hTmp1 = plot([1;1] * trig.SF{blk}(:,1)' / H.sampleRate, [0 freqLim],'-g','LineWidth',3);
            hTmp2 = plot([1;1] * trig.SB{blk}(:,1)' / H.sampleRate, [0 freqLim],'-r','LineWidth',3);
            h = [hTmp1(1) hTmp2(1)];
            lH = legend(h,{'Extenstion','Retraction'});
            set(lH)
            colorbar
            % set(gca,'ydir','normal','ylim',trft_frequencies([1 end]),'xlim',trftOfBlock_time([1 end])/1000)
            set(gca,'ydir','normal','ylim',trft_frequencies([1 freqLimInd]),'xlim',trftOfBlock_time([1 end])/H.sampleRate)
            xlabel('Time (s)')
            ylabel('Frequency (Hz)')
            title([H.datasetName, ' Block ' num2str(useBlocks(blk)), ' ', eeggroup], 'FontSize', 36)
        end
    end
    %%

    [SF_trft.data,~,trialTime] = trigeredAmplitude(EEG.data, ...
                                                             trft_baseMean, ...
                                                             trig.SF, ...
                                                             fftTrialPast, ...
                                                             fftTrialFuture, ...
                                                             fftWinSize, ...
                                                             fftStep, ...
                                                             H.sampleRate, ...
                                                             fftWinFunction, ...
                                                             [], ...
                                                             true);

    SF_trft.mean = zeros(size(SF_trft.data,1),size(SF_trft.data,2),H.noOfEl);
    SF_trft.std = zeros(size(SF_trft.data,1),size(SF_trft.data,2),H.noOfEl);
    SF_trft.noTrials = size(SF_trft.data,4);
    for ch = 1:H.noOfEl
        SF_trft.mean(:,:,ch) = nanmean(SF_trft.data(:,:,ch,:),4);
        SF_trft.std(:,:,ch) = nanstd(SF_trft.data(:,:,ch,:),0,4);
    end

    [SB_trft.data,~,trialTime] = trigeredAmplitude(EEG.data, ...
                                                             trft_baseMean, ...
                                                             trig.SB, ...
                                                             fftTrialPast, ...
                                                             fftTrialFuture, ...
                                                             fftWinSize, ...
                                                             fftStep, ...
                                                             H.sampleRate, ...
                                                             fftWinFunction, ...
                                                             [], ...
                                                             true);

    SB_trft.mean = zeros(size(SB_trft.data,1),size(SB_trft.data,2),H.noOfEl);
    SB_trft.std = zeros(size(SB_trft.data,1),size(SB_trft.data,2),H.noOfEl);
    SB_trft.noTrials = size(SB_trft.data,4);
    for ch = 1:H.noOfEl
        SB_trft.mean(:,:,ch) = nanmean(SB_trft.data(:,:,ch,:),4);
        SB_trft.std(:,:,ch) = nanstd(SB_trft.data(:,:,ch,:),0,4);
    end

    %% Calculate trft mean by electrodes
    eeggroups = {'frontal', 'left', 'center', 'right', 'parietal'};
    Neegs = length(eeggroups);
    SF_trft.mean_eeggroups = zeros(size(SF_trft.mean,1), size(SF_trft.mean, 2), Neegs);
    SB_trft.mean_eeggroups = zeros(size(SB_trft.mean,1), size(SB_trft.mean, 2), Neegs);
    for i=1:Neegs
        eeggroup = eeggroups{i};
        chaninds = getEEGChannelInds(X(1).chanlocs, eeggroup);
        SF_trft.mean_eegroups(:,:,i) = squeeze(mean(SF_trft.mean(:,:,chaninds), 3));
        SB_trft.mean_eegroups(:,:,i) = squeeze(mean(SB_trft.mean(:,:,chaninds), 3));
    end


    %%
    maxF = 25;


    clear param;
    param.yAxis = trft_frequencies(1:maxF);
    param.xAxis = trialTime/1000;
    param.xlabel = 'Time relative to the cue (s)';
    param.ylabel = 'Frequency (Hz)';
    param.plotTriggerMarker = true;
    % param.triggersAt = [0 0.3];

    %param.clabel = 'Normalized spectral amplitude (a.u.)';
    param.useColorbar = 'true';
    param.design = {[], 1, [];
                    2,  3,  4;
                    [], 5, [];};
    param.titleFont = 35;
    titleString = 'Reach Extension - average spectral amplitude normalized to mean';
    param.pictureTitle = titleString;
    %figHandle = plotMultiarray(SF_trft.mean_eegroups(1:maxF,:,:),param);


    titleString = 'Reach retraction - average spectral amplitude normalized to mean';
    param.pictureTitle = titleString;
    %figHandle = plotMultiarray(SB_trft.mean_eegroups(1:maxF,:,:),param);

    %%

    trft_baseNorm = baselineAmplitude_SE(EEG.data, ...
                                        trft_baseMean, ...
                                        [], ...   % TRIGGER.baseBeta_SE, ...
                                        baseOffset, ...
                                        fftWinSize, ...
                                        fftBaseStep, ...
                                        H.sampleRate, ...
                                        fftWinFunction);


    %% Looking for optimal bands

    maxFreq = 25;

    trft_bandSnr_SF = zeros(maxFreq,maxFreq,length(trialTime),H.noOfEl);
    trft_bandSnr_SB = zeros(maxFreq,maxFreq,length(trialTime),H.noOfEl);
    for b1 = 1:maxFreq
        cueBand_SF = 0;
        cueBand_SB = 0;
        baseBand = 0;
        for b2 = b1:maxFreq
            cueBand_SF = (b2-b1)/(b2-b1+1)*cueBand_SF + squeeze(SF_trft.data(b2,:,:,:))/(b2-b1+1);
            cueBand_SB = (b2-b1)/(b2-b1+1)*cueBand_SB + squeeze(SB_trft.data(b2,:,:,:))/(b2-b1+1);
            baseBand = (b2-b1)/(b2-b1+1)*baseBand + squeeze(trft_baseNorm(b2,:,:))/(b2-b1+1);

            trft_bandSnr_SF(b1,b2,:,:) = abs(mean(cueBand_SF,3)-repmat(mean(baseBand,2)',[size(cueBand_SF,1) 1]))./...
                                         (std(cueBand_SF,0,3)+repmat(std(baseBand,0,2)',[size(cueBand_SF,1) 1]));

            trft_bandSnr_SB(b1,b2,:,:) = abs(mean(cueBand_SB,3)-repmat(mean(baseBand,2)',[size(cueBand_SB,1) 1]))./...
                                         (std(cueBand_SB,0,3)+repmat(std(baseBand,0,2)',[size(cueBand_SB,1) 1]));

            disp([num2str(b1) '/' num2str(maxFreq) '; '...
                  num2str(b2) '/' num2str(maxFreq)]);
        end
    end


    %% Plotting SNR for all bands

    clear param;
    param.yAxis = trft_frequencies(1:maxFreq);
    param.xAxis = trft_frequencies(1:maxFreq);
    param.xlabel = 'Band Top Frequency (Hz)';
    param.ylabel = 'Band Bottom Frequency (Hz)';
    param.clabel = {'SNR'};
    param.useColorbar = 'true';
    param.design = H.electrodeConfig;


    param.pictureTitle = 'Reach extension: frequency band SNR';
    max_trft_bandSnr_SF = squeeze(max(trft_bandSnr_SF,[],3));
    % plotMultiarray(max_trft_bandSnr_SF,param)

    param.pictureTitle = 'Reach retraction: frequency band SNR';
    max_trft_bandSnr_SB = squeeze(max(trft_bandSnr_SB,[],3));
    % plotMultiarray(max_trft_bandSnr_SB, param)

    %% Calculate average SNR across electrodes
    eeggroups = {'frontal', 'left', 'center', 'right', 'parietal'};
    Neegs = length(eeggroups);
    trft_bandSnr_SF_eeggroups = zeros(maxFreq,maxFreq,Neegs);
    trft_bandSnr_SB_eeggroups=  zeros(maxFreq,maxFreq,Neegs);
    for i=1:Neegs
        eeggroup = eeggroups{i};
        chaninds = getEEGChannelInds(X(1).chanlocs, eeggroup);
        trft_bandSnr_SF_eeggroups(:,:,i) = squeeze(mean(max_trft_bandSnr_SF(:,:,chaninds), 3));
        trft_bandSnr_SB_eeggroups(:,:,i) = squeeze(mean(max_trft_bandSnr_SB(:,:,chaninds), 3));
    end

    param.design = {[], 1, [];
                    2,  3,  4;
                    [], 5, [];};

    param.pictureTitle = sprintf('Subject %d, reach extension: frequency band SNR', subject);
    h = plotMultiarray(trft_bandSnr_SF_eeggroups,param);
    fname = sprintf('SNR_extension_subject%d', subject);
    saveas(h, fname, 'png');

    param.pictureTitle = sprintf('Subject %d, reach retraction: frequency band SNR', subject);
    h = plotMultiarray(trft_bandSnr_SB_eeggroups,param);
    fname = sprintf('SNR_retraction_subject%d', subject);
    saveas(h, fname, 'png');

    % %% Calculate SNR
    % for blk = 1:H.noOfBlocks
    %     [trftOfBlock,~,trftOfBlock_time] = baselineAmplitude_SE(EEG.data, ...
    %         trft_baseMean, ... % [], ...
    %         [], ... %TRIGGER.baseStamp_SE_noReal, ...
    %         baseOffset, ...
    %         fftWinSize, ...
    %         pltStep, ...
    %         H.sampleRate, ...
    %         fftWinFunction, ...
    %         blk);
    %     
    %     [SNR] = Calculate_SNR(trftOfBlock, ...
    %         trftOfBlock_time/H.sampleRate, ...
    %         [1] * trig.SF{blk}' / H.sampleRate, ...
    %         fftTrialPast/H.sampleRate, ...
    %         fftTrialFuture/H.sampleRate);
    %     
    %     SNR_blk(:,:,blk) = SNR;
    %     
    % end
    % 
    % SNR_mean = mean(SNR_blk,3);
    % 
    % 
    % 

    %% Looking for optimal bands using SQNA

    SF_trft.data = sqrt(SF_trft.data);
    SB_trft.data = sqrt(SB_trft.data);
    trft_baseNorm = sqrt(trft_baseNorm);

    maxFreq = 25;

    trft_bandSnr_SF = zeros(maxFreq,maxFreq,length(trialTime),H.noOfEl);
    trft_bandSnr_SB = zeros(maxFreq,maxFreq,length(trialTime),H.noOfEl);
    for b1 = 1:maxFreq
        cueBand_SF = 0;
        cueBand_SB = 0;
        baseBand = 0;
        for b2 = b1:maxFreq
            cueBand_SF = (b2-b1)/(b2-b1+1)*cueBand_SF + squeeze(SF_trft.data(b2,:,:,:))/(b2-b1+1);
            cueBand_SB = (b2-b1)/(b2-b1+1)*cueBand_SB + squeeze(SB_trft.data(b2,:,:,:))/(b2-b1+1);
            baseBand = (b2-b1)/(b2-b1+1)*baseBand + squeeze(trft_baseNorm(b2,:,:))/(b2-b1+1);

            trft_bandSnr_SF(b1,b2,:,:) = abs(mean(cueBand_SF,3)-repmat(mean(baseBand,2)',[size(cueBand_SF,1) 1]))./...
                                         (std(cueBand_SF,0,3)+repmat(std(baseBand,0,2)',[size(cueBand_SF,1) 1]));

            trft_bandSnr_SB(b1,b2,:,:) = abs(mean(cueBand_SB,3)-repmat(mean(baseBand,2)',[size(cueBand_SB,1) 1]))./...
                                         (std(cueBand_SB,0,3)+repmat(std(baseBand,0,2)',[size(cueBand_SB,1) 1]));

            disp([num2str(b1) '/' num2str(maxFreq) '; '...
                  num2str(b2) '/' num2str(maxFreq)]);
        end
    end


    %% Plotting SNR for all bands

    clear param;
    param.yAxis = trft_frequencies(1:maxFreq);
    param.xAxis = trft_frequencies(1:maxFreq);
    param.xlabel = 'Band Top Frequency (Hz)';
    param.ylabel = 'Band Bottom Frequency (Hz)';
    param.clabel = {'SNR'};
    param.useColorbar = 'true';
    param.design = H.electrodeConfig;


    param.pictureTitle = 'Reach extension: frequency band SNR';
    max_trft_bandSnr_SF = squeeze(max(trft_bandSnr_SF,[],3));
    % plotMultiarray(max_trft_bandSnr_SF,param)

    param.pictureTitle = 'Reach retraction: frequency band SNR';
    max_trft_bandSnr_SB = squeeze(max(trft_bandSnr_SB,[],3));
    % plotMultiarray(max_trft_bandSnr_SB, param)

    %% Calculate average SNR across electrodes
    eeggroups = {'frontal', 'left', 'center', 'right', 'parietal'};
    Neegs = length(eeggroups);
    trft_bandSnr_SF_eeggroups = zeros(maxFreq,maxFreq,Neegs);
    trft_bandSnr_SB_eeggroups=  zeros(maxFreq,maxFreq,Neegs);
    for i=1:Neegs
        eeggroup = eeggroups{i};
        chaninds = getEEGChannelInds(X(1).chanlocs, eeggroup);
        trft_bandSnr_SF_eeggroups(:,:,i) = squeeze(mean(max_trft_bandSnr_SF(:,:,chaninds), 3));
        trft_bandSnr_SB_eeggroups(:,:,i) = squeeze(mean(max_trft_bandSnr_SB(:,:,chaninds), 3));
    end

    param.design = {[], 1, [];
                    2,  3,  4;
                    [], 5, [];};

    param.pictureTitle = sprintf('Subject %d, reach extension: frequency band SNR SQNA', subject);
    h = plotMultiarray(trft_bandSnr_SF_eeggroups,param);
    fname = sprintf('SNR_extension_SQNA_subject%d', subject);
    saveas(h, fname, 'png');

    param.pictureTitle = sprintf('Subject %d, reach retraction: frequency band SNR SQNA', subject);
    h = plotMultiarray(trft_bandSnr_SB_eeggroups,param);
    fname = sprintf('SNR_retraction_SQNA_subject%d', subject);
    saveas(h, fname, 'png');
end




% 
