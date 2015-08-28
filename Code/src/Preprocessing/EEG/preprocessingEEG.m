function preprocessingEEG(subject, exptype)
%%PREPROCESSINGEEG - performs preprocessing on EEG data from experiments
%  preprocessingEEG()
%
%  Write a full length description
%
%  PARAMETERS:
%    exptype: string - specifies the type of experiment.  Options are:
%            'Healthy' 'IAF', 'LT', 'MJ', and 'OAF'
%
%  RETURNS:
%    specify return variables if any
%
%  @ June 2015, Sean Bittner   sbittner@andrew.cmu.edu
%
%    - adapted from file preprocessing_reaching_passive_HM_NEW 
%      written by Aurelie Sadaka-Stephan

%% Set parameters
if isnumeric(subject)
    subject = num2str(subject);
end
channels = {'Fp1', 'Fpz', 'Fp2', 'F1', 'Fz', 'F2', ...                          % Frontal
            'F5', 'F3', 'FC5', 'FC3', 'C5', 'C3', 'CP5', 'CP3', 'P5', 'P3', ... % Left
            'F4', 'F6', 'FC4', 'FC6', 'C4', 'C6', 'CP4', 'CP6', 'P4', 'P6', ... % Right
            'FC1', 'FCz', 'FC2', 'C1', 'Cz', 'C2', 'CP1', 'CPz', 'CP2', ...     % Central
            'P1', 'Pz', 'P2', 'POz', ...                                        % Parietal
            'F7', 'F8', 'T7', 'T8', 'P7', 'P8', 'O1', 'O2', 'Oz'};
if (strcmp(exptype, 'Healthy'))
    expstr = 'HM';
else
    expstr = exptype;
end
block_events = cell(1,10);
block_names = cell(1,10);
for i = 1:12
    block_events{i} = sprintf('%s_%d_event', expstr, i);
    if (strcmp(subject, '3') || strcmp(subject, '4'))
        block_names{i} = sprintf('%s%d.bdf', exptype, i);
    else
        block_names{i} = sprintf('%s_%d.bdf', exptype, i);
    end
end
sr=2048;
sr_new=1000;
low_cf=1;
high_cf=50;
filt_order=8;
rule_ICA=20;
prestimulus=0.2;
n_blocks=length(block_events);
    
%% Setup file I/O
if (ispc) % Setting up file I/O for windows
    % Set filepath locations
    baseDir = ['E:', filesep, 'Sean', filesep];
    dataDir = [baseDir, 'Data', filesep];
    eegPreprocessDir = [baseDir, 'Code', filesep, 'src', filesep, ...
                       'Preprocessing', filesep, 'EEG', filesep 'Subject', ...
                        subject, filesep];
    eegDir = [dataDir, 'EEG', filesep];
    eventDir = [dataDir, 'Event', filesep];
    %addpath(genpath(baseDir)); %Add code from all subfolders to matlab path
    
    ch_file=[eegDir, filesep, 'channels_location(64).ced']; % to change 

    eventSubjectDir = [eventDir, 'Subject_' subject, filesep];
    eegSubjectDir = [eegDir, 'Subject', subject, filesep];

    namefolder=[eegSubjectDir 'Reaching HM'];
    if (exist(namefolder) ~= 7)
        mkdir(namefolder);
    end
end


%% Preprocess EEG data
% Determine which indices must be kept according to specified electrodes.
chlocs=readlocs(ch_file);
nchannels = length(channels);
channels_to_keep = zeros(1,nchannels);
j = 1;
for i=1:64
    if (nnz(ismember(channels, chlocs(i).labels)))
        channels_to_keep(j) = i;
        j = j + 1;
    end
end
chlocs = chlocs(channels_to_keep);
        
    

% Data loading
eeglab;

% Extra channels removing and channel locations reading
fname_removed = [eegPreprocessDir, exptype, '_all_ch_removed1.mat'];
if (exist(fname_removed, 'file') == 2)
    fprintf('Already have selected data\n');
    load(fname_removed);
else
    all_original=struct([]);
    for i=1:n_blocks
        disp(['Loading data... ' block_names{i}]);
        block_fname = [eegSubjectDir, block_names{i}];
        data=pop_biosig(block_fname);
        data.setname=['Sub' subject '_' block_names{i} '_original'];
        all_original=struct([all_original data]);
    end
    all_ch_removed1=struct([]);
    for i=1:n_blocks
        disp(['Only keeping channels specified by 5 regions ... ' block_names{i}]);
        data=all_original(i);
        data.data=data.data(channels_to_keep,:);
        data.nbchan=nchannels;
        data.chanlocs=chlocs;
        data.setname=['Sub' subject '_' block_names{i} '_ch_removed1'];
        all_ch_removed1=struct([all_ch_removed1 data]);
    end
    save(fname_removed,'all_ch_removed1','-v7.3');
end

%% Data filtering
fname_filtered = [eegPreprocessDir, exptype, '_all_filtered.mat'];
if (exist(fname_filtered, 'file') == 2)
    fprintf('Already have filtered data\n');
    load(fname_filtered);
else
    all_filtered=struct([]);
    h1=fdesign.highpass('N,F3dB',filt_order,low_cf,sr);
    d1=design(h1,'Butter');
    h2=fdesign.lowpass('N,F3dB',filt_order,high_cf,sr);
    d2=design(h2,'Butter');
%     freqz(d1) % don't need to see these filters
%     freqz(d2)
    for i=1:n_blocks
        disp(['Filtering data... '  block_names{i}]);
        data=all_ch_removed1(i);
        data.data = double(data.data);
        data.data=data.data-repmat(mean(data.data,2),[1 size(data.data,2)]);
        for nn=1:size(data.data,1)
            data.data(nn,:)=filtfilt(d1.sosMatrix,d1.ScaleValues,...
                data.data(nn,:));
            data.data(nn,:)=filtfilt(d2.sosMatrix,d2.ScaleValues,...
                data.data(nn,:));
        end
        data.setname=['Sub' subject '_' block_names{i} '_filtered'];
        all_filtered=struct([all_filtered data]);
    end
    save(fname_filtered,'all_filtered','-v7.3');
end

%% Data resampling
all_resampled=struct([]);
fname_resampled = [eegPreprocessDir, exptype, '_all_resampled.mat'];
if (exist(fname_resampled, 'file') == 2)
    fprintf('Already have resampled data\n');
    load(fname_resampled);
else
    for i=1:n_blocks        
        disp(['Resampling data... ' block_names{i}]);
        %data=pop_resample(all_filtered(i),sr_new);     Pop resample is the spawn of matlab satan.
        block = all_filtered(i);
        data_resampled = resample(block.data', sr_new, sr)';
        time_resampled = resample(block.times', sr_new, sr)';
        block.data = double(data_resampled);
        block.srate = sr_new;
        block.pnts = size(data_resampled, 2);
        block.times = double(time_resampled);
        data.setname=['Sub' subject '_' block_names{i} '_resampled'];       
        all_resampled=struct([all_resampled block]);   
    end
    save(fname_resampled,'all_resampled','-v7.3');
end

%% Channels visual removal
all_ch_removed2=struct([]);
fname_badchanremoved = [eegPreprocessDir, exptype, '_all_ch_removed2.mat'];
if (exist(fname_badchanremoved, 'file') == 2)
    fprintf('Already have data with bad channels removed\n');
    load(fname_badchanremoved);
else
    for ii=1:n_blocks
        disp(['Channels visual removal... ' block_names{ii}]);
        data=all_resampled(ii);
        figure
        eegplot(data.data,'srate',sr_new);
    end
    bad_chans=input('Bad channels (vector): ');
    for ii=1:n_blocks
        data=all_resampled(ii);
        good_chans=setdiff(1:nchannels,bad_chans);
        n_channels=length(good_chans);

        data.data=data.data(good_chans,:);
        data.nbchan=n_channels;
        data.chanlocs=data.chanlocs(good_chans);
        data.setname=['Sub' subject '_' block_names{ii} '_ch_removed2'];
        all_ch_removed2=struct([all_ch_removed2 data]);
    end
    save(fname_badchanremoved,'all_ch_removed2');
end

%% Common average re-referencing
fname_reref = [eegPreprocessDir, exptype, '_all_rereferenced.mat'];
if (exist(fname_reref, 'file') == 2)
    fprintf('Already have re-referenced data\n');
    load(fname_reref);
else
    all_rereferenced=struct([]);
    for i=1:n_blocks
        disp(['Re-referencing data... ' block_names{i}]);
        data=all_ch_removed2(i);
        data.data=data.data-repmat(mean(data.data,1),[size(data.data,1) 1]);
        data.setname=['Sub' subject '_' block_names{i} '_rereferenced'];
        all_rereferenced=struct([all_rereferenced data]);
    end
    
    % After rereferencing, the data are now linearly dependent.
    % I will remove the channel Oz so we can perform ICA to remove
    % blinking components from the data   
    chlocs = all_rereferenced(1).chanlocs;
    nchannels = length(chlocs);
    for i=1:nchannels
        if strcmp(chlocs(i).labels, 'Oz')
            ind = i;
        end
    end
    channels_to_keep = [1:(ind -1), (ind + 1):nchannels];
    for i=1:n_blocks
        data=all_rereferenced(i);
        data.data=data.data(channels_to_keep,:);
        data.nbchan=nchannels-1;
        data.chanlocs=chlocs(channels_to_keep);
        all_rereferenced(i) = data;
    end
    save(fname_reref,'all_rereferenced');
end

%% Independet Component Analisys (ICA): Infomax
all_ICA1=struct([]);
n_comp=zeros(1,n_blocks);

fname_ICA1 = [eegPreprocessDir, exptype, '_all_ICA1.mat'];
if (exist(fname_ICA1, 'file') == 2)
    fprintf('Already have ICA data\n');
    load(fname_ICA1);
else
    for i=1:n_blocks
        disp(['Running Infomax ICA... ' block_names{i}]);
        block=all_rereferenced(i);
        tmpdata=block.data-repmat(mean(block.data,2),[1 size(block.data,2)]);
        n_comp(i)=min([rank(tmpdata),floor(sqrt((size(block.data,2)/...
            rule_ICA)))]);
        data=pop_runica(block,'icatype','runica','extended',1,'pca',...
            n_comp(i),'stop',1e-7);
        data.setname=['Sub' subject '_' block_names{i} '_ICA1'];
        all_ICA1=struct([all_ICA1 data]);
        % Condition number revision
    end
    save(fname_ICA1,'all_ICA1', 'n_comp');
    for i=1:n_blocks
        cond_num=cond(all_ICA1(i).icaweights);
        disp([block_names{i} ' -> ' num2str(cond_num)]);
    end
end



%% Remove ICA components via visual inspection
all_IC_visual_removal1=struct([]);
fname_ICA_remove = [eegPreprocessDir, exptype, '_all_IC_visual_removal1.mat'];
if (exist(fname_ICA_remove) == 2)
    fprintf('Already have ICA component removed data\n');
    load(fname_ICA_remove);
else
    for i=1:n_blocks
        disp(['Independent Components Inspection... ' block_names{i}]);
        data=all_ICA1(i);
        pop_topoplot(data,0,1:n_comp(i),'Independent Components');
        eegplot(eeg_getdatact(data,'component',1:n_comp(i),'reshape','2d'));
        bad_IC=input('Bad componets (vector): ');
        if(~isempty(bad_IC))
            projection=data.icawinv(:,bad_IC)*eeg_getdatact(data,...
                'component',bad_IC,'reshape','2d');
            data.data=data.data-projection;
            good_IC=setdiff(1:size(data.icaweights,1),bad_IC);
            data.icawinv=data.icawinv(:,good_IC);
            data.icaweights=data.icaweights(good_IC,:);
        end
        data.reject=[];
        data.setname=['Sub' subject '_' block_names{i} '_IC_visual_removal1'];
        all_IC_visual_removal1=struct([all_IC_visual_removal1 data]);
    end
    save(fname_ICA_remove, 'all_IC_visual_removal1');
end
Trials = all_IC_visual_removal1;
resfname = [eegSubjectDir, sprintf('Preprocessed%s.mat', exptype)];
save(resfname, 'Trials');

end

