function Y = extractEEGSignalFast(data, Event, srate_new, targets, extraction_type, varargin)
%EXTRACTEEGSIGNAL - Get preprocessed eeg signal for a particular reach
%   extractEEGSignalFast(data, Event, srate_new, targets, extraction_type, varargin)
%   This function is useful for extracting the epoch from eeg data, when
%   the eeg data for a particular experiment have already been preloaded,
%   along with the corresponding data structure.  This function is useful
%   when multiple epochs need to be extracted from the same EEG file.
%      
%   If you choose the 'uniform length' option for the extraction_type
%   parameter, you will have to enter two additional arguements.  The first
%   is the pre move onset buffer in seconds, and the second is the post
%   move onset buffer in seconds.  In the following example, a pre move
%   buffer of 0.7 seconds and a post move buffer of 2.51 seconds are
%   specified.
%   Ex:
%   extractEEGSignalFast(...,'uniform length', .70, 2.51);
%
%  PARAMETERS
%    data - samples x nchannels matrix of eeg data
%    Event - event data structure corresponding to data matrix
%
%    channels - string specifying particular EEG channel or EEG group
%      Ex: 'Fp1', 'Cz', etc.                                    (single)
%      Ex: 'parietal', 'fronatl', 'center', 'left', or 'right'  (group)
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
%  RETURNS
%    Y - an m x n matrix where m is the number of electrodes and n is the
%        number of samples. Each row is a vector of preprocessed EEG
%        corresponding to the specified reaching epoch.
%
%  @ 2015 Sean Bittner   sbittner@andrew.cmu.edu
%

srate = 2048;

switch extraction_type
    case 'reach duration'
        [begin_ind, end_ind] = getReachIndsReachDuration(Event.EEG, Event, targets, ...
                                                     srate, srate_new);
    case 'uniform length'
        if (length(varargin) > 1)
            premove_buf = varargin{1};
            postmove_buf = varargin{2};
        else
            fprintf('Error: Must supply premove buffer and postmove buffer arguments.\n');
        end
        [begin_ind, end_ind] = getReachIndsUniformLength(Event.EEG, Event, targets, ...
                                        srate, srate_new, premove_buf, postmove_buf);
    case 'percentiles'
        if (length(varargin) > 0)
            windowsize = varargin{1};
        else
            fprintf('Error: Must supply window size argument.\n');
        end
        [begin_ind, end_ind] = getReachIndsPercentiles(Event.EEG, Event, targets, srate, srate_new, windowsize);
end
Y = data(begin_ind:end_ind,:);

end

