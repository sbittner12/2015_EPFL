function [channel_ind] = getEMGChannelInd(EMG, musclestr)
%GETEMGCHANNELIND - Provides data matrix channel index from string spec.
%   Based on a muscle string specification the function returns the
%   corresponding index of the data matrix associated with the channel
%   location.
%
%  PARAMETERS
%    EMG - EMG struct of corresponding data
%
%    musclestr - string specifying particular EMG channel
%      Ex: 'DANT', 'BICS', etc.
% 
%  RETURNS
%    channel_ind - integer index to extract from data matrix
%
%  @ 2015 Sean Bittner   sbittner@andrew.cmu.edu
%
    % Determine muscle index
    emgDir = ['E:',filesep,'Sean', filesep, 'Data', filesep, 'EMG', filesep];
    load([emgDir, 'EMG_labels.mat']);
    nchanlocs = length(correct_labels);
    if (strcmp(musclestr, 'all'))
        channel_ind = 1:nchanlocs;
    else
        for i=1:nchanlocs
            if (strcmp(musclestr, correct_labels{i}))
                channel_ind = i;
                break;
            elseif (i == nchanlocs)
                fprintf('Error: could not find muscle %s.\n', musclestr);
                channel_ind = -1
            end
        end
    end


end


