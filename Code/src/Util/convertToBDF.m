function [] = convertToBDF(subject, trial, target)
%%CONVERTTOBDF - converts EEG single epoch .mat files to .bdf files
%   convertToBDF(subject, trial, target)
%
%   PARAMETERS:
%     subject - scalar : subject id
%     trial   - scalar : trial id
%     target  - scalar : target id
%
%   RETURNS:
%     none
%
%   @ 2015 Sean Bittner sbittner@andrew.cmu.edu   
%
    baseDir = ['E:', filesep, 'Sean', filesep, 'Data', filesep];
    eegDir = [baseDir, 'EEG', filesep];

    eegSubjectEpochMatDir = [eegDir, sprintf('Subject%d', subject), ...
                          filesep, 'SingleEpochMat', filesep];
    eegSubjectEpochBdfDir = [eegDir, sprintf('Subject%d', subject), ...
                          filesep, 'SingleEpochBdf', filesep];
    if (exist(eegSubjectEpochBdfDir, 'dir') ~= 7)
        mkdir(eegSubjectEpochBdfDir);
    end              
    eegfnameprefix = sprintf('OAF_%d_%d', trial, target);

    load([eegSubjectEpochMatDir, sprintf('%s.mat', eegfnameprefix)]);
    resfname = [eegSubjectEpochBdfDir, sprintf('%s.bdf', eegfnameprefix)];
    pop_writeeeg(EEG, resfname, 'TYPE', 'BDF');

end

