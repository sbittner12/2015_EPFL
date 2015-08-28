%run_3way_mbPLS_subjects
clear all

subjects = [3:6,8:9];
plotavg = 1;
exptype = 'OAF';
channels = 'all';
nsubjects = length(subjects);
% choose decomposition method
eegMode = 'eegSpectrogram'; %{eegSpectrogram, eegCoherence}
emgMode = 'emgRaw'; %{imfAmp, emgAmp)

correct_labels = {'TRAPS', 'TRAPM', 'DANT', 'DMED', 'DPOS', 'PEC', 'INFRA',...
                  'LAT', 'RHO', 'BICL', 'BICS' 'TRILAT', 'TRILONG', 'BRAC', 'PRO'};

% choose decomposition method
eegMode = 'eegSpectrogram'; %{eegSpectrogram, eegCoherence}
emgMode = 'emgAmp'; %{imfAmp, emgAmp)

baseDir = ['E:', filesep, 'Sean', filesep];
dataDir = [baseDir, 'Data', filesep];
resDir = [baseDir, 'Results', filesep, 'PartialLeastSquares', filesep];
eegDir = [dataDir, 'EEG', filesep];
emgDir = [dataDir, 'EMG', filesep];

resfname = [resDir, sprintf('%s_%s_AllSubjects.mat', eegMode, emgMode)];
resampledspectfname = [resDir, 'Spectrogram_AllSubjects.mat'];
resampledemgfname = [resDir, 'EMG_AllSubjects.mat'];


Xdata_groups = cell(1, nsubjects);
Ydata_groups = cell(1, nsubjects);
% Load the Spectrograms and Raw EMG
for i=1:nsubjects
    subject = subjects(i);
    fprintf('Subject %d\n', subject);
    fprintf('Loading EEG and EMG signals.\n')

    eegSubjectDir = [eegDir, sprintf('Subject%d', subject), filesep];
    spectfname = [eegSubjectDir, 'Spectrogram_AllTrials.mat'];
    
    fprintf('Loading subject %d spectrograms ...', subject);
    X = load(spectfname);
    Xdata_groups{i} = X.Xdata_group;
    fprintf('done.\n');
    
    fprintf('Loading subject %d emg ...', subject);
    emgSubjectDir = [emgDir, sprintf('Subject%d', subject), filesep];
    emgfname = [emgSubjectDir, 'EMG_AllTrials.mat'];
    Y = load(emgfname);
    Ydata_groups{i} = Y.Ydata_group;
    fprintf('done.\n');
end

Xdata_group = cell(1,nsubjects);
Ydata_group = cell(1,nsubjects);
for i=1:nsubjects
    i
    Xdata_group{i}
    Xpls = trialBlockstoSubjectBlock(Xdata_groups{i});
    Xdata_group{i} = Xpls;
    Ypls = trialBlockstoSubjectBlock(Ydata_groups{i});
    Ydata_group{i} = Ypls;
    size(Xpls)
    size(Ypls)
end

    % perform MBPLS decomposition
numFactor = 4;
fprintf('running mbpls');
[Tb,W1b,W2b,Wt,Tt,Ub,Qb,Wu,Tu, ssx, ssy, W1b_temp, W2b_temp, Tu_orig] = ...
    MBtriPLS2_cmmnW1b(Xdata_group, Ydata_group, numFactor, size(Xdata_group{1}, 1)*1e-12);

EEG = loadEEG(3, 'OAF', 'all', 1000);
F = (3:.5:50)';


for f=1:numFactor
    figure; 

    subplot(3, 4, [1,2]);
    plot(F, abs(W1b{1}(:, f))); title('Subject 3 spectral signature');
    xlim([min(F) max(F)]); 
    xlabel('Frequency (Hz)'); 

    subplot(3, 4, [5,6]);
    nsamples = size(Tt,1);
    plot(1:nsamples, unitVariance( Tt(:, f).'), 'r', ... 
         1:nsamples, unitVariance(Tu(:, f).'), 'b'); %, ...
         %1:nsamples, unitVariance(Tu_orig(:,f).'), 'g'); 
    legend('EEG', 'EMG');%, 'original EMG activation');
    title('Temporal Signature');
    xlim([1 nsamples]); 
    
    if (plotavg)
        Qbavg = averageQb(Qb, f);
        subplot(3, 4, [9,10]), bar(Qbavg), title('Average muscle synergy');
        ax = gca;
        set(ax, 'XTickLabel', correct_labels, 'FontSize', 6);
    else
        [Qb_dotprods, Qb_colors] = rankQb(Qb, f);
        [~, ind] = max(Qb_dotprods);
        subplot(3, 4, [9,10]), bar(Qb{ind}(:,f), 'FaceColor', Qb_colors{ind}), title('Best fit synergy');
        ax = gca;
        set(ax, 'XTickLabel', correct_labels, 'FontSize', 6);
        subplot(3, 4, [13,14]), bar(Qb_dotprods), title('Dot product with best fit synergy');
        ax = gca;
        set(ax, 'XTickLabel', {'3', '4', '5', '6', '7', '8', '9'}, 'FontSize', 12);
        xlabel('subject');
    end
    
    for i=1:nsubjects
        ind = 4*ceil(i/2) - mod(i,2);
        subplot(3, 4, ind);
        subject = subjects(i);
        topoplot(W2b{i}(:,f), EEG(1).chanlocs);
        title(sprintf('Subject %d', subject));
    end
        
    suptitle(sprintf('Component %i', f));

    
%         for b = 1:ntrials
end
