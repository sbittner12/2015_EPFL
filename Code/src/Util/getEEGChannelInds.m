function [channel_inds] = getEEGChannelInds(chanlocs, channels)
%GETEEGCHANNELINDS - Provides data matrix channel index from string spec.
%   Based on a string specification (which can be that of a single
%   electrode or a group of electrodes), the function returns the
%   corresponding indices of the data matrix associated with the channel
%   locations.
%
%  PARAMETERS
%    chanlocs - property of an eeglab .bdf file which specifies the eeg
%               electrode layout of the corresponding data matrix.
%    channels - string specifying particular EEG channel or EEG group
%      Ex: 'Fp1', 'Cz', etc.                                    (single)
%      Ex: 'parietal', 'fronatl', 'center', 'left', or 'right'  (group) 
% 
%  RETURNS
%    channel_inds - list of integers of indices to extract from data matrix
%
%  @ 2015 Sean Bittner   sbittner@andrew.cmu.edu
%
    channelGroupsDir = ['E:',filesep,'Sean', filesep, 'Code', filesep, ...
                        'src', filesep, 'Coherence', ...
                        filesep, 'ChannelGroupings', filesep];
                    
    nchanlocs = length(chanlocs);
    % Find which channels are to be extracted based on channels string
    switch channels
        case 'frontal'
            F = load([channelGroupsDir, 'frontal_Fang.mat']);
            chanlocs_to_get = F.chanlocs;
        case 'left'
            F = load([channelGroupsDir, 'left_Fang.mat']);
            chanlocs_to_get = F.chanlocs;
        case 'right'
            F = load([channelGroupsDir, 'right_Fang.mat']);
            chanlocs_to_get = F.chanlocs;
        case 'center'
            F = load([channelGroupsDir, 'center_Fang.mat']);
            chanlocs_to_get = F.chanlocs;
        case 'parietal'
            F = load([channelGroupsDir, 'parietal_Fang.mat']);
            chanlocs_to_get = F.chanlocs;
        case 'all'
            F = load([channelGroupsDir, 'all.mat']);
            chanlocs_to_get = F.chanlocs;
        case 'all_Fang'
            F = load([channelGroupsDir, 'all_Fang.mat']);
            chanlocs_to_get = F.chanlocs;
        otherwise % assume that the channel string is one channel to extract
            chanlocs_to_get = {channels};
    end

    % Find which indices of the data matrix should be kept based on chanlocs
    nchanlocs_to_get = length(chanlocs_to_get);
    ind = 1;
    for i=1:nchanlocs_to_get
        chanloc_to_get = chanlocs_to_get{i};
        for j=1:nchanlocs
            if (strcmp(chanloc_to_get, chanlocs(j).labels))
                %fprintf('Using channel location %s.\n',chanloc_to_get);
                channel_inds(ind) = j;
                ind = ind + 1;
                break;
            elseif (j == nchanlocs)
              %fprintf('Channel location %s not found in preprocessed data.\n',chanloc_to_get);
            end
        end
    end     


end

