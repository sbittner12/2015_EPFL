function Y = extractEMGSignalFast(emgdata, Event, targets, extraction_type, varargin)
%EXTRACTEMGSIGNAL - Get preprocessed EMG signal for a particular reach
%   extractEMGSignal(subject, exptype, musclestr, trial, target, extraction_type, varargin)
%   This function is useful for extracting epochs of EMG data that
%   correspond to particular reaches and performing preprocessing on the 
%   extracted signal. By specifying the subject and experiment type you can
%   extract specified EMG electrode data.  You can selct one target, or 
%   choose to show the entire trial with the 'all' string.
%   Ex:
%   extractEMGSignal(3, 'OAF', 'DANT', 1, 1, 'reach duration');   
%      
%   If you choose the 'uniform length' option for the extraction_type
%   parameter, you will have to enter two additional arguements.  The first
%   is the pre move onset buffer in seconds, and the second is the post
%   move onset buffer in seconds.  In the following example, a pre move
%   buffer of 0.7 seconds and a post move buffer of 2.51 seconds are
%   specified.
%   Ex:
%   extractEMGSignal(3,'OAF','DANT',1,1,'uniform length', .70, 2.51);
%
%  PARAMETERS
%    subject - integer or number string which identifies the subject
%      Ex: '3', or 3
%    exptype - string for experiment type based on exoskeleton function
%      Ex: 'Healthy', 'IAF', 'LT', 'MJ', or 'OAF'
%    musclestr - string specifying particular EMG channel
%      Ex: 'DANT', 'BICS', etc.
%    trial - integer specifying trial number
%
%    target - integer specifying target or string 'all'
%              Ex: 1, or 'all' for the full trial signal
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
%             Note: only one target can be extracted in this case.
%     
%    varargin - when the extraction_type is 'uniform length', the first
%             parameter must be the specified pre movement onset buffer in 
%             seconds while the second parameter must be the post movement
%             onset buffer in seconds
%
%  RETURNS
%    Y - a 1 x n vector of preprocessed EMG data corresponding to the 
%        specified reaching epoch where n is the number of samples.
%
%  @ 2015 Sean Bittner   sbittner@andrew.cmu.edu
%
srate_new = 1000;
srate = 3000;


switch extraction_type
   case 'reach duration'
        [begin_ind, end_ind] = getReachIndsReachDuration(Event.EMG, Event, targets, srate, srate_new);
    case 'uniform length'
        if (length(varargin) > 1)
            premove_buf = varargin{1};
            postmove_buf = varargin{2};
        else
            fprintf('Error: Must supply premove buffer and postmove buffer arguments.\n');
        end
        [begin_ind, end_ind] = getReachIndsUniformLength(Event.EMG, Event, targets, ...
                                     srate, srate_new, premove_buf, postmove_buf);
    case 'percentiles'
        if (length(varargin) > 0)
            windowsize = varargin{1};
        else
            fprintf('Error: Must supply window size argument.\n');
        end
        [begin_ind, end_ind] = getReachIndsPercentiles(Event.EMG, Event, targets, srate, srate_new, windowsize);
end

Y = emgdata(begin_ind:end_ind,:);

end

