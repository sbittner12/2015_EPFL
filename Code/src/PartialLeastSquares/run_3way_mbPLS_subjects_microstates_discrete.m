%run_3way_mbPLS_subjects
clear all

subjects = [3,4,5,6,7,8,9];
exptype = 'OAF';
trials = 1:12;
targets = 1:8;
channels = 'all';
Neeg = 46;
Nemg = 15;
plotavg = 1;
nsubjects = length(subjects);
ntrials = length(trials);
ntargets = length(targets);
nmicrostates = 6;
green = [0,153,0]/255;
% choose decomposition method
emgMode = 'emgRaw'; %{imfAmp, emgAmp)

correct_labels = {'TRAPS', 'TRAPM', 'DANT', 'DMED', 'DPOS', 'PEC', 'INFRA',...
                  'LAT', 'RHO', 'BICL', 'BICS' 'TRILAT', 'TRILONG', 'BRAC', 'PRO'};

min_epoch =  getShortestEpoch();
nsamples = min_epoch*ntrials*ntargets;

% choose decomposition method
eegMode = 'eegSpectrogram'; %{eegSpectrogram, eegCoherence}
emgMode = 'emgAmp'; %{imfAmp, emgAmp)

baseDir = ['E:', filesep, 'Sean', filesep];
dataDir = [baseDir, 'Data', filesep];
resDir = [baseDir, 'Results', filesep, 'PLS', filesep];
eegDir = [dataDir, 'EEG', filesep];
emgDir = [dataDir, 'EMG', filesep];
resfname = [eegDir, 'MicrostateDiscrete_AllTrials.mat'];


if (exist(resfname, 'file') == 2)
    load(resfname)
else
    eegSubjectData = zeros(min_epoch*ntrials*ntargets, Neeg);
    emgSubjectData = zeros(min_epoch*ntrials*ntargets, Nemg);
    % Load the eeg and emg data
    Xrawdata_group = cell(nsubjects,1);
    Ydata_group = cell(nsubjects,1);
    targetinds = zeros(12*min_epoch,8);
    for i=1:nsubjects
        subject = subjects(i);
        eegSubjectDir = [eegDir, sprintf('Subject%d', subject), filesep];
        emgSubjectDir = [emgDir, sprintf('Subject%d', subject), filesep];
        for j=1:ntrials
            trial = trials(j);
            %fprintf('trial %d\n', j);

            for k=1:ntargets
                target = targets(k);
                %fprintf('target %d\n', k);
                ind1 = min_epoch*((j-1)*ntargets+k-1) + 1;
                ind2 = min_epoch*((j-1)*ntargets+k);

                %update target inds
                tind1 = min_epoch*(j-1)+1;
                tind2 = min_epoch*(j);
                targetinds(tind1:tind2, k) = (ind1:ind2)';


                eegfname = [eegSubjectDir, 'SingleEpochMat', filesep, sprintf('%s_%d_%d.mat', exptype, trial, target)];
                load(eegfname);
                eegData = EEG.data;
                eegSubjectData(ind1:ind2, :) = eegData';

                emgepochfname = [emgSubjectDir, 'SingleEpochMat', filesep, sprintf('%s_%d_%d.mat', exptype, trial, target)];
                load(emgepochfname);
                emgData = EMG.data;
                emgSubjectData(ind1:ind2, :) = emgData';
            end
        end
    %     Xpls = nprocess(eegSubjectData, [1, 0], [0, 1], [], [], 1, 1);
    %     Ypls = nprocess(emgSubjectData, [1, 0], [0, 1], [], [], 1, 1);

        Xrawdata_group{i} = eegSubjectData;
        Ydata_group{i} = nprocess(emgSubjectData, [1, 0], [0, 1], [], [], 1, 1);
    end

     
    % Load the microstates
    fprintf('Loading microstates...\n');
    microstates = loadMicrostates(subjects, nmicrostates);
    [microstates_clustered, ref_maps, groups] = clusterMicrostates(microstates, subjects,1);
    
    % Fit the microstates to the input EEG
    fprintf('Fitting microstates to the data EEG using correlation coefficients...\n');
    Xdata_group = cell(nsubjects, 1);
    for i=1:nsubjects
        subject = subjects(i);
        fprintf('Subject %d\n', subject);
        subjectEEG = Xrawdata_group{i};
        subjectMicrostates = microstates_clustered{i};
        microstate_activations = zeros(nsamples, nmicrostates);
        for j=1:nmicrostates
            fprintf('microstate %d\n', j);
            for k=1:nsamples
                microstate_activations(k,j) = corr2(subjectEEG(k,:), subjectMicrostates(:,j)');
            end
        end
        % Discretize the microstate activations
        microstates_discrete = discretizeMicrostates(microstate_activations, nsamples);
        Xpls = nprocess(microstates_discrete, [1,0], [0,1], [], [], 1, 1);
        Xdata_group{i} = Xpls;
    end
    
    save(resfname, 'Xdata_group', 'Ydata_group', 'microstates', 'microstates_clustered', 'ref_maps', 'groups', 'targetinds')
end

EEG = loadEEG(3, 'OAF', 'all', 2048);
chanlocs = EEG(1).chanlocs;
plotReferenceMicrostates(microstates, ref_maps, chanlocs);

%% perform MBPLS decomposition
numFactor = 4;
fprintf('running mbpls\n');
[Tb,Pb,Wb, Wb_reproj, Wt,Tt,Ub,Qb,Wu,Tu, ssx, ssy, Wridge] = ...
    MBbiPLS2(Xdata_group, Ydata_group, numFactor, size(Xdata_group{1}, 1)*1e-12);


%% Plot the components
plottype = 2;

switch plottype
    case 1
        rows = 4+nsubjects;
    case 2
        rows = 5;
end

for f=1:numFactor
    figure; 
    Tt_target1 = crossTrialAverage(Tt(targetinds(:,1),f), min_epoch, ntrials);
    Tu_target1 = crossTrialAverage(Tu(targetinds(:,1),f), min_epoch, ntrials);
    Tt_target3 = crossTrialAverage(Tt(targetinds(:,3),f), min_epoch, ntrials);
    Tu_target3 = crossTrialAverage(Tu(targetinds(:,3),f), min_epoch, ntrials);
    Tt_target5 = crossTrialAverage(Tt(targetinds(:,5),f), min_epoch, ntrials);
    Tu_target5 = crossTrialAverage(Tu(targetinds(:,5),f), min_epoch, ntrials);
    Tt_target7 = crossTrialAverage(Tt(targetinds(:,7),f), min_epoch, ntrials);
    Tu_target7 = crossTrialAverage(Tu(targetinds(:,7),f), min_epoch, ntrials);
    
    
    subplot(rows, 4, [2,3]);
    samples = 1:min_epoch;
    plot(samples, unitVariance(Tt_target1), 'r', samples, unitVariance(Tu_target1), 'b'); 
    title('Temporal signature for target 1');
    xlim([1 min_epoch]); 
    
    subplot(rows, 4, [7,8]);
    plot(samples, unitVariance(Tt_target3), 'r', samples, unitVariance(Tu_target3), 'b'); 
    title('Temporal signature for target 3');
    xlim([1 min_epoch]); 
    
    subplot(rows, 4, [10,11]);
    plot(samples, unitVariance(Tt_target5), 'r', samples, unitVariance(Tu_target5), 'b'); 
    title('Temporal signature for target 5');
    xlim([1 min_epoch]); 
    
    subplot(rows, 4, [5,6]);
    plot(samples, unitVariance(Tt_target7), 'r', samples, unitVariance(Tu_target7), 'b'); 
    title('Temporal signature for target 7');
    xlim([1 min_epoch]); 
    
    % Plot best fit microstate
    [Pb_dotprods, Pb_colors] = rankQb(Pb, f);
    [~, ind] = max(Pb_dotprods);
    subplot(rows, 4, [13,14]), bar(Pb{ind}(:,f), 'FaceColor', green), title('Best fit microstate comb');
    xlabel('microstates');
    
    % Plot best fit muscle synergy
    [Qb_dotprods, Qb_colors] = rankQb(Qb, f);
    if (plotavg)
        Qbavg = averageQb(Qb, f);
        subplot(rows, 4, [15,16]), bar(Qbavg, 'FaceColor', green), title('Average muscle synergy');
        ax = gca;
        set(ax, 'XTickLabel', correct_labels, 'FontSize', 6);
    else
        [~, ind] = max(Qb_dotprods);
        subplot(rows, 4, [15,16]), bar(Qb{ind}(:,f), 'FaceColor', green), title('Best fit synergy');
        ax = gca;
        set(ax, 'XTickLabel', correct_labels, 'FontSize', 6);
    end
    

    switch plottype
        case 1
            % Plot microstate combination and muscle synergy for each group
            for i=1:nsubjects
                subject = subjects(i);
                ind = 16+4*i;

                subplot(rows, 4, [ind-3, ind-2]);
                bar(Pb{i}(:,f), 'FaceColor', Pb_colors{i}), title(sprintf('Subject %d microstate comb', subject));
                xlabel('microstates');

                subplot(rows, 4, [ind-1, ind]);
                bar(Qb{i}(:,f), 'FaceColor', Qb_colors{i}), title(sprintf('Subject %d muscle synergy', subject));
            end
        case 2
            subplot(rows,4, [17,18]), bar(Pb_dotprods), title('Dot product with best fit microstate comb');
            xlabel('subjects'); 
            
            subplot(rows,4, [19,20]), bar(Qb_dotprods), title('Dot product with best fit synergy');
            xlabel('subjects');
    end
        
    suptitle(sprintf('Component %i', f));
end
