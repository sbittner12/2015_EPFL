% function [FS_trft,FO_trft] = Analysis_rat_608_SNR()

clear all;

% Directory that holds the data from the session
useBlocks = [1:4 8 10:11];

plotDiagnostics = true;

H.noOfBlocks = length(useBlocks);
H.datasetName = 'Rat 610';
H.defectiveEl = 2;
H.noOfEl = 8;
H.CARch = [];
H.electrodeConfig = {8, 3, []; ...
                     7, 5, 1; ...
                     6, 4, 2};

electrodeList = 1:H.noOfEl;
H.sampleRate = 2.441406250000000e+04;

%% Collecting raw data file names
% These files contain 30kHz broadband recordings (ns5)
ECoG.data = cell(length(useBlocks),1);
trig.move = cell(length(useBlocks),1);
trig.calm = cell(length(useBlocks),1);
trig.task = cell(length(useBlocks),1);
for blk = 1:length(useBlocks)
%     load(['E:\ECoG_Version1407\608\P21\608_P21_RW' num2str(useBlocks(blk)) '_Neu_noREF_cut.mat']);
%     load(['E:\ECoG_Version1407\608\P7\608_P7_RW' num2str(useBlocks(blk)) '_Neu_noREF_cut.mat']);
    load(['E:\Data\DataECoG\Analysis\610\P55\610_P55_RW_' num2str(useBlocks(blk)) '_Neu_noREF_cut.mat']);
%     load(['E:\ECoG_Version1407\608\P55_LD\608_P55_LD' num2str(useBlocks(blk)) '_Neu_noREF_cut.mat']);

    %ECoG.data{blk} = data_cut(:,1:8) - repmat(data_cut(:,9),[1 8]);
    ECoG.data{blk} = data_cut(:,1:8);              
%     load(['E:\ECoG_Version1407\608\P21\608_P21_RW' num2str(useBlocks(blk)) '_Neu_noREF_trig.mat']);
%     load(['E:\ECoG_Version1407\608\P7\608_P7_RW' num2str(useBlocks(blk)) '_Neu_noREF_trig.mat']);
    load(['E:\Data\DataECoG\Analysis\610\P55\610_P55_RW_' num2str(useBlocks(blk)) '_Neu_noREF_trig.mat']);
%     load(['E:\ECoG_Version1407\608\P55_LD\608_P55_LD' num2str(useBlocks(blk)) '_Neu_noREF_trig.mat']);

%     trig.FO{blk} = trigger_TDT.FO - trigger_TDT.FO(1);
%     trig.FS{blk} = trigger_TDT.FS - trigger_TDT.FO(1);
    
    trig.FO{blk} = trigger_TDT.FO;
    trig.FS{blk} = trigger_TDT.FS;

end

%% Doing the STFT analysis

trialPast = round(0.1 * H.sampleRate);
trialFuture = round(0.3 * H.sampleRate);
trialTime = (-trialPast:trialFuture)/H.sampleRate;
fftWinSize = ceil(H.sampleRate/10);
fftWinFunction = hamming(fftWinSize);
fftStep = floor(H.sampleRate/40);
fftBaseStep = floor(fftWinSize/2);
fftTrialPast = trialPast + floor(fftWinSize/2);
fftTrialFuture = trialFuture + floor(fftWinSize/2);
offsetBefore = fftTrialPast;
offsetAfter = fftTrialFuture;

usedElectrodes = 1:H.noOfEl;



%% Plotting the trft of the whole block together with the cues


freqLim = 1000;
freqLimInd = find(trft_frequencies < freqLim,1,'last');

pltStep = floor(fftWinSize/16);
% meanEl = setdiff(1:H.noOfEl,H.noisyEl);

meanEl = [1 7];

for blk = 1:H.noOfBlocks
    [trftOfBlock,~,trftOfBlock_time] = baselineAmplitude_SE(ECoG.data, ...
                                                        trft_baseMean, ... % [], ...  
                                                        [], ... %TRIGGER.baseStamp_SE_noReal, ...
                                                        baseOffset, ...
                                                        fftWinSize, ...
                                                        pltStep, ...
                                                        H.sampleRate, ...
                                                        fftWinFunction, ...
                                                        blk);

    trftOfBlockMean = squeeze(mean(trftOfBlock(:,meanEl,:),2));

    figure
    hold on
    imagesc(trftOfBlock_time/H.sampleRate, trft_frequencies,trftOfBlockMean)
    hTmp1 = plot([1;1] * trig.FO{blk}' / H.sampleRate, [0 freqLim],'--w','LineWidth',1);
    hTmp2 = plot([1;1] * trig.FS{blk}' / H.sampleRate, [0 freqLim],'-w','LineWidth',1);
    h = [hTmp1(1) hTmp2(1)];
    lH = legend(h,{'FO','FS'});
    set(lH, 'Color', 'red')
    colorbar
    % set(gca,'ydir','normal','ylim',trft_frequencies([1 end]),'xlim',trftOfBlock_time([1 end])/1000)
    set(gca,'ydir','normal','ylim',trft_frequencies([1 freqLimInd]),'xlim',trftOfBlock_time([1 end])/H.sampleRate)
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title([H.datasetName ' block ' num2str(useBlocks(blk)) ' - STFT amplitude normalized to baseline, averaged over channels'])
end

%%

[FS_trft.data,~,trialTime] = trigeredAmplitude(ECoG.data, ...
                                                         trft_baseMean, ...
                                                         trig.FS, ...
                                                         fftTrialPast, ...
                                                         fftTrialFuture, ...
                                                         fftWinSize, ...
                                                         fftStep, ...
                                                         H.sampleRate, ...
                                                         fftWinFunction, ...
                                                         [], ...
                                                         true);

FS_trft.mean = zeros(size(FS_trft.data,1),size(FS_trft.data,2),H.noOfEl);
FS_trft.std = zeros(size(FS_trft.data,1),size(FS_trft.data,2),H.noOfEl);
FS_trft.noTrials = size(FS_trft.data,4);
for ch = 1:H.noOfEl
    FS_trft.mean(:,:,ch) = nanmean(FS_trft.data(:,:,ch,:),4);
    FS_trft.std(:,:,ch) = nanstd(FS_trft.data(:,:,ch,:),0,4);
end

[FO_trft.data,~,trialTime] = trigeredAmplitude(ECoG.data, ...
                                                         trft_baseMean, ...
                                                         trig.FO, ...
                                                         fftTrialPast, ...
                                                         fftTrialFuture, ...
                                                         fftWinSize, ...
                                                         fftStep, ...
                                                         H.sampleRate, ...
                                                         fftWinFunction, ...
                                                         [], ...
                                                         true);

FO_trft.mean = zeros(size(FO_trft.data,1),size(FO_trft.data,2),H.noOfEl);
FO_trft.std = zeros(size(FO_trft.data,1),size(FO_trft.data,2),H.noOfEl);
FO_trft.noTrials = size(FO_trft.data,4);
for ch = 1:H.noOfEl
    FO_trft.mean(:,:,ch) = nanmean(FO_trft.data(:,:,ch,:),4);
    FO_trft.std(:,:,ch) = nanstd(FO_trft.data(:,:,ch,:),0,4);
end

%%

datasetName = '607';

maxF = 60;

clear param;
param.yAxis = trft_frequencies(1:maxF);
param.xAxis = trialTime/1000;
param.xlabel = 'Time relative to the cue (s)';
param.ylabel = 'Frequency (Hz)';
param.plotTriggerMarker = true;
% param.triggersAt = [0 0.3];

param.clabel = 'Normalized spectral amplitude (a.u.)';
param.useColorbar = 'true';
param.design = H.electrodeConfig;
titleString = ', Foot strike, LFPs response to stimulation - average spectral amplitude normalized to mean';
param.pictureTitle = [datasetName titleString];
figHandle = plotMultiarray(FS_trft.mean(1:maxF,:,:),param);



titleString = ', Foot off, LFPs response to stimulation - average spectral amplitude normalized to mean';
param.pictureTitle = [datasetName titleString];
figHandle = plotMultiarray(FO_trft.mean(1:maxF,:,:),param);

%% Single trials

band = 19:22;

clear param;
param.xAxis = trialTime/1000;
param.xlabel = 'Time relative to the cue (s)';
param.ylabel = 'Trials';
param.clabel = 'Spike rate (a.u.)';
param.useColorbar = 'true';
param.design = H.electrodeConfig;

param.yAxis = 1:size(FS_trft.data,4);
titleString = ', single trial spike rate responses to online detected Foot strike events followed by stimulation';
param.pictureTitle = [datasetName titleString];
figHandle = plotMultiarray(permute(squeeze(mean(FS_trft.data(band,:,:,:),1)),[3 1 2]),param);

param.yAxis = 1:size(FO_trft.data,4);
titleString = ', single trial spike rate responses to online detected Foot strike events followed by stimulation';
param.pictureTitle = [datasetName titleString];
figHandle = plotMultiarray(permute(squeeze(mean(FO_trft.data(band,:,:,:),1)),[3 1 2]),param);


%%

trft_baseNorm = baselineAmplitude_SE(ECoG.data, ...
                                    trft_baseMean, ...
                                    [], ...   % TRIGGER.baseBeta_SE, ...
                                    baseOffset, ...
                                    fftWinSize, ...
                                    fftBaseStep, ...
                                    H.sampleRate, ...
                                    fftWinFunction);
                                             

%% Looking for optimal bands

maxFreq = 102;

trft_bandSnr_FS = zeros(maxFreq,maxFreq,length(trialTime),H.noOfEl);
trft_bandSnr_FO = zeros(maxFreq,maxFreq,length(trialTime),H.noOfEl);
for b1 = 1:maxFreq
    cueBand_FS = 0;
    cueBand_FO = 0;
    baseBand = 0;
    for b2 = b1:maxFreq
        cueBand_FS = (b2-b1)/(b2-b1+1)*cueBand_FS + squeeze(FS_trft.data(b2,:,:,:))/(b2-b1+1);
        cueBand_FO = (b2-b1)/(b2-b1+1)*cueBand_FO + squeeze(FO_trft.data(b2,:,:,:))/(b2-b1+1);
        baseBand = (b2-b1)/(b2-b1+1)*baseBand + squeeze(trft_baseNorm(b2,:,:))/(b2-b1+1);

        trft_bandSnr_FS(b1,b2,:,:) = abs(mean(cueBand_FS,3)-repmat(mean(baseBand,2)',[size(cueBand_FS,1) 1]))./...
                                     (std(cueBand_FS,0,3)+repmat(std(baseBand,0,2)',[size(cueBand_FS,1) 1]));

        trft_bandSnr_FO(b1,b2,:,:) = abs(mean(cueBand_FO,3)-repmat(mean(baseBand,2)',[size(cueBand_FO,1) 1]))./...
                                     (std(cueBand_FO,0,3)+repmat(std(baseBand,0,2)',[size(cueBand_FO,1) 1]));

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


param.pictureTitle = 'Foot Strike: frequency band SNR';
plotMultiarray(squeeze(max(trft_bandSnr_FS,[],3)),param)

param.pictureTitle = 'Foot Off: frequency band SNR';
plotMultiarray(squeeze(max(trft_bandSnr_FO,[],3)),param)