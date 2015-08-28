%%CONVERTALLTOBDF - converts all EEG single epoch files to .bdf files
%
%   By specifying the subjects, trials, and targets of the epochs, convert
%   the .mat files to .bdf files which can be used for EEG microstates
%   extraction.
%
%   @ 2015 Sean Bittner sbittner@andrew.cmu.edu
%

subjects = 3:9
trials = 1:12;
targets = 1:8;

for subject = subjects
    fprintf('subject %d\n', subject);
    for trial = trials
        fprintf('trial %d\n', trial);
        for target = targets
            convertToBDF(subject, trial, target);
        end
    end
end