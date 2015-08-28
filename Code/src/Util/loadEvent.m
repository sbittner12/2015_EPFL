function [event] = loadEvent(subject, exptype, trial)
%LOADEVENT Summary of this function goes here
%   loadEvent(subject, exptype, trial)
%   Detailed explanation goes here

% Setup file I/O
if isnumeric(subject)
    subject = num2str(subject);
end
dataDir = ['E:',filesep,'Sean', filesep, 'Data', filesep];
subjectEventDir = [dataDir, 'Event', filesep, 'Subject_', subject,filesep]; 

eventstr = getEventstr(exptype);

% Find indices to extract proper EEG signal from data
if (trial < 1 || 12 < trial)
    fprintf('Error: trial must be integer 1 through 12\n');
end

eventfname = [subjectEventDir, sprintf('%s_%d_event.mat', eventstr, trial)];
if (exist(eventfname, 'file') == 2)
    load(eventfname);
    event = Event;
else
    fprintf('Event for subject %s, exptype %s, trial %d does not exist.\n', ...
             subject, exptype, trial);
    event = {};
end

end

