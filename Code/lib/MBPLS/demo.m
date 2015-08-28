% demo
clear;

addpath(genpath(pwd));
dummyData = load('demo_data'); % load dummy EEG and EMG data

numSubj = 5;

Xdata_group = cell(1, numSubj);
Ydata_group = cell(1, numSubj);

% choose decomposition method
eegMode = 'rawEEG'; %{rawEEG, eegCoherence}
emgMode = 'emgAmp'; %{imfAmp, emgAmp)


for i = 1:numSubj
    
    % load i-th subject's EEG and EMG data
    % here we create dummy data from eegData and emgData variables 
    eegData = dummyData.eegData + randn(size(dummyData.eegData));
    emgData = dummyData.emgData + randn(size(dummyData.emgData));
    
    % remove channel mean and set to unit variance
    eegData = unitVariance(removeMean(eegData));
    emgData = unitVariance(removeMean(emgData));
    
    switch emgMode % extract temporal features of EMG data
        case 'imfAmp'  % extract amplitudes of EMD components
            imfData = [];
            for ch = 1:size(emgData, 1)
                imf = emd2(emgData(ch, :), 'maxmodes', 5, 'display', 0);
                imfData = [imfData; imf(1:3, :)];
            end

            imfAmp = zeros(size(imfData));
            for ch = 1:size(imfData)
                imfAmp(ch, :) = est_amp( imfData(ch, :), 2, 60 );
            end
            emgComp = imfAmp;
            
        case 'emgAmp' % extract amplitude of raw EMG data
            emgAmp = zeros(size(emgData));
            for ch = 1:size(emgData)
                emgAmp(ch, :) = est_amp( emgData(ch, :), 2, 60 );
            end      
            emgComp = emgAmp;
    end   
    
    switch eegMode % extract temporal features of EEG data
        case 'rawEEG' % spectrogram of raw EEG data
            eegSData = [];
            eeg_sample_rate = 250; 
            WL = eeg_sample_rate*1.5;
            for ch = 1:size(eegData)
                [S,F,T] = spectrogram(eegData(ch, :), WL, ...
                    WL-5, 3:0.5:50, eeg_sample_rate, 'yaxis');

                eegSData(:, :, ch) = (abs((S).'));
            end
            
            Xpls=nprocess(eegSData,[1 0 0], [0 0 0], [], [], 1, 1);

        case 'eegCoherence' % coherence between EEG channels
            WL = 500;
            T = [];
            curpos = 0; ll = 0;
            while curpos+WL < size(eegData, 2)
                T(ll+1)= curpos/eeg_sample_rate;
                curpos = curpos + shift_len;
                ll = ll +1; 
            end
            eegCoh = NaN(ll, 100, nchoosek(size(eegData, 1), 2));
            counter = 1; 
            for m = 1:size(eegData, 1)
                for n = m+1:size(eegData, 1)   
                    for z = 1:ll
                        index = (z-1)*shift_len + 1 :  ((z-1)*shift_len + WL);
                        [Cxy,F] = mscohere(eegData(m, index),eegData(n, index), ...
                            hamming(100), 100-shift_len, ...
                            eeg_sample_rate*2,eeg_sample_rate); 
                        eegCoh(z, :, counter) = Cxy((F>0)& (F<=50));
                        F = F((F>0)& (F<=50));
                    end
                    counter = counter+1;
                end
            end
    
            Xpls=nprocess(eegCoh,[1 0 0], [0 0 0], [], [], 1, 1);
    end

    % temporally align eeg and emg data
    size(emgComp)
    size(Xpls)
    Ypls=nprocess( resample(emgComp.', length(Xpls), length(emgComp)),...
        [1 0], [0 1], [], [], 1, 1);
    
    Xdata_group{i} = Xpls;
    Ydata_group{i} = Ypls;

end


% perform MBPLS decomposition
numFactor = 10;
[Tb,W1b,W2b,Wt,Tt,Ub,Qb,Wu,Tu, ssx, ssy, W1b_temp, W2b_temp] = ...
    MBtriPLS2_cmmnW1b(Xdata_group,Ydata_group, numFactor,  size(Xdata_group{1}, 1)*1e-12);

% display temporal and spectral signatures of each component
for f = 1:numFactor
    figure; 
    
    subplot(2, 2, [1 2]);
    plot(F, abs(W1b{1}(:, f))); title('Spectral Signature');
    xlim([min(F) max(F)]); 
    xlabel('Frequency (Hz)'); 
    
    subplot(2, 2, [3 4]);
    plot(T, unitVariance( Tt(:, f).'), 'r', ... 
         T, unitVariance(Tu(:, f).'), 'b'); 
    legend('EEG', 'EMG'); title('Temporal Signature');
    xlim([min(T) max(T)]); 
    xlabel('Time (sec)'); 
    

    suptitle(sprintf('Component %i', f));
end
