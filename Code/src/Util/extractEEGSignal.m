function Y = extractEEGSignal(subject, exptype, channels, trial, tragets, extraction_type, varargin)
%EXTRACTEEGSIGNAL - Get preprocessed eeg signal for a particular reach
%   extractEEGSignal(subject, exptype, channels, trial, targets, extraction_type, varargin)
%   This function is useful for finding epochs of preprocessed EEG data
%   that correspond to particular reaches.  By specifying the subject and 
%   experiment type you can extract specified EEG electrode data.  You can
%   also extract groups of electrodes using group descriptors like
%   'frontal' or 'parietal'.  You can choose to select multiple targets,
%   which will return the smallest window of data including all targets.
%   Ex:
%   extractEEGSignal(3, 'OAF', 'frontal', 1, 1, 'reach duration');   
%      
%   If you choose the 'uniform length' option for the extraction_type
%   parameter, you will have to enter two additional arguements.  The first
%   is the pre move onset buffer in seconds, and the second is the post
%   move onset buffer in seconds.  In the following example, a pre move
%   buffer of 0.7 seconds and a post move buffer of 2.51 seconds are
%   specified.
%   Ex:
%   extractEEGSignal(3,'OAF','frontal',1,1,'uniform length', .70, 2.51);
%
%  PARAMETERS
%    subject - integer or number string which identifies the subject
%      Ex: '3', or 3
%    exptype - string for experiment type based on exoskeleton function
%      Ex: 'Healthy', 'IAF', 'LT', 'MJ', or 'OAF'
%    channels - string specifying particular EEG channel or EEG group
%      Ex: 'Fp1', 'Cz', etc.                                    (single)
%      Ex: 'parietal', 'fronatl', 'center', 'left', or 'right'  (group)
%    trial - integer specifying trial number
%
%    targets - list of integers specifying targets to show
%              Ex: 1, or 1:8 for the full trial signal
%
%    extraction_type - string specifying which samples to extract for epoch
%      Ex: 'reach duration' - takes all samples for center -> target and
%             from target -> center.  Specifically this is from the trigger
%             index to the second corresponding stop index
%      Ex: 'uniform length' - takes samples before and after motion onset
%             specified by two varargin parameters. The first parameter
%             specifies the length in seconds prior to motion onset to
%             extract samples, and the second parameter specifies the
%             seconds after motion onset to keep samples.  Motion onset is
%             determined by the first start signal of a particular epoch.
%             Typically the pre movement and post movement onset buffers
%             should be calculated using the getMovementBuffers function in
%             the util directory.
%      Ex: 'percentiles' - takes the 0th, 10th, 20th, ... 100th percentiles
%      of the start-end period.  Requires extra parameter windowsize
%     
%    varargin - when the extraction_type is 'uniform length', the first
%             parameter must be the specified pre movement onset buffer in 
%             seconds while the second parameter must be the post movement
%             onset buffer in seconds.  When the extraction_type is
%             'percentiles', the first additional parameter must be the
%             windowsize
%
%
%  RETURNS
%    Y - an m x n matrix where m is the number of electrodes and n is the
%        number of samples. Each row is a vector of preprocessed EEG
%        corresponding to the specified reaching epoch.
%
%  @ 2015 Sean Bittner   sbittner@andrew.cmu.edu
%

srate = 2048;

% Setup file I/O
if isnumeric(subject)
    subject = num2str(subject);
end
dataDir = ['E:',filesep,'Sean', filesep, 'Data', filesep];
subjectEEGDir = [dataDir, 'EEG', filesep, 'Subject', subject, filesep];
subjectEventDir = [dataDir, 'Event', filesep, 'Subject_', subject,filesep]; 

eventstr = getEventstr(exptype);
                  
fname =  [subjectEEGDir, 'Preprocessed', exptype, '.mat'];
if (~(exist(fname, 'file') == 2))
    fprintf('Error: file %s does not exist.\n', fname);
    return;
else
    load(fname);
end

% Determine channel indices
chanlocs = Trials(trial).chanlocs;
channel_inds = getEEGChannelInds(chanlocs, channels);

    
% Find indices to extract proper EEG signal from data
if (trial < 1 || 12 < trial)
    fprintf('Error: trial must be integer 1 through 12\n');
end

srate_new = Trials(trial).srate;
data = Trials(trial).data;

eventfname = [subjectEventDir, sprintf('%s_%d_event', eventstr, trial)];
load(eventfname);

if (strcmp(extraction_type, 'uniform length'))
    if (length(varargin) > 1)
        premove_buf = varargin{1};
        postmove_buf = varargin{2};
    end
    [begin_ind, end_ind] = getReachIndsUniformLength(Event.EEG, targets, ...
                                    srate, srate_new, premove_buf, postmove_buf);
elseif (strcmp(extraction_type, 'reach duration'))
    [begin_target_ind, end_target_ind] = parseTargets(targets, Event);
    [begin_ind, end_ind] = getReachIndsReachDuration(Event.EEG, begin_target_ind, end_target_ind, ...
                                                     srate, srate_new);
end
Y = data(channel_inds, begin_ind:end_ind);

end

