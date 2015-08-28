function [pls] = trialBlockstoSubjectBlock(data_group)
%TRIALBLOCKSTOSUBJECTBLOCK Summary of this function goes here
%   Detailed explanation goes here
ndatagroups = length(data_group);
pls = [];
for i=1:ndatagroups
    pls = cat(1, pls, data_group{i});
end

end

