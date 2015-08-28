function [chanlocs] = loadChanlocs(channels)
%LOADCHANLOCS Summary of this function goes here
%   Detailed explanation goes here
EEG = loadEEG(3,'OAF',channels, 2048);
chanlocs = EEG(1).chanlocs;

end

