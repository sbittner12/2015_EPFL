function [N] = getNumEEGs(channels)
%GETNUMEEGS Summary of this function goes here
%   Detailed explanation goes here
channelGroupsDir = ['E:',filesep,'Sean', filesep, 'Code', filesep, ...
                    'src', filesep, 'Coherence', ...
                    filesep, 'ChannelGroupings', filesep];
                    
switch channels
    case 'frontal'
        F = load([channelGroupsDir, 'frontal_Fang.mat']);
    case 'left'
        F = load([channelGroupsDir, 'left_Fang.mat']);
    case 'right'
        F = load([channelGroupsDir, 'right_Fang.mat']);
    case 'center'
        F = load([channelGroupsDir, 'center_Fang.mat']);
    case 'parietal'
        F = load([channelGroupsDir, 'parietal_Fang.mat']);
    case 'all'
        F = load([channelGroupsDir, 'all.mat']);
    otherwise % assume that the channel string is one channel to extract
        fprintf('Error: EEG channel group %s not recognized.\n', channels);
        N = -1;
        return;
end
N = length(F.chanlocs);

end

