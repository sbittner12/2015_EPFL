function [microstates] = loadMicrostates(subjects, nmicrostates)
%LOADMICROSTATES Summary of this function goes here
%   Detailed explanation goes here

if (nmicrostates < 0)
    nmicrostates = 1;
elseif (nmicrostates > 12)
    nmicrostates = 12;
end

nsubjects = length(subjects);

baseDir = ['E:', filesep, 'Sean', filesep];
dataDir = [baseDir, 'Data', filesep];
eegDir = [dataDir, 'EEG', filesep];

microstates = cell(nsubjects, 1);
for i=1:nsubjects
    subject = subjects(i);
    microstatefname = [eegDir, sprintf('Subject%d', subject), filesep, 'SingleEpochBdf', ...
                       filesep, 'Seg OAF.OAF_1_1', filesep, ...
                       sprintf('Seg OAF.OAF_1_1.0%d.(0%d).ep', nmicrostates, nmicrostates)];
    microstates{i} = load(microstatefname).';
end

