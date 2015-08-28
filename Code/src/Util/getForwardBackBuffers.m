function [preforward_buf, postforward_buf, prebackward_buf, postbackward_buf] = getForwardBackBuffers(subjects, varargin)
%SUBJECTINFO Summary of this function goes here
%   Detailed explanation goes here
% Setup file I/O
EEG_Fs = 2048;
nsubjects = length(subjects);
ntrials = 12;

nvarargin = length(varargin);
if (nvarargin > 0)
    exptypes = {varargin{1}};
    nexptypes = 1;
else
    exptypes = {'HM', 'IAF', 'LT', 'MJ', 'OAF'};
    nexptypes = length(exptypes);
end

dataDir = ['E:',filesep,'Sean', filesep, 'Data', filesep];

preforward_buf = 10000*ones(nexptypes,1);
postforward_buf = 10000*ones(nexptypes,1);
prebackward_buf = 10000*ones(nexptypes,1);
postbackward_buf = 10000*ones(nexptypes,1);

for k = 1:nsubjects
    subject = subjects(k);
    subjectEventDir = [dataDir, 'Event', filesep, sprintf('Subject_%s', subject), filesep];
    
    for i=1:nexptypes
        exptype = exptypes{i};
        fprintf('******** %s ********\n', exptype);
        min_preforward = [];
        min_postforward = [];
        min_prebackward = [];
        min_postbackward = [];
        ind = 1;
        for j=1:ntrials
            eventfname = [subjectEventDir, sprintf('%s_%d_event.mat',exptype, j)];
            if (exist(eventfname, 'file') == 2)
                load(eventfname);
                triggers = parseTriggers(Event.EEG.Trigger, Event.EEG.End(end));
                % Find min pre forward movement period
                preforward = (Event.EEG.Start(1:2:end-1) - triggers);
                min_preforward(ind) = min(preforward);
                % Find min post forward movement period
                postforward = (Event.EEG.End(1:2:end-1) - Event.EEG.Start(1:2:end-1));
                min_postforward(ind) = min(postforward);
                % Find min pre backward movement period
                prebackward = (Event.EEG.Start(2:2:end) - Event.EEG.End(1:2:end-1));
                min_prebackward(ind) = min(prebackward);
                % Find min post backward movement period
                postbackward = (Event.EEG.End(2:2:end) - Event.EEG.Start(2:2:end));
                min_postbackward(ind) = min(postbackward);
                ind = ind + 1;
            else
                display('file dne');
            end
        end
        if (preforward_buf(i) > min(min_preforward))
            preforward_buf(i) = min(min_preforward);
        end
        if (postforward_buf(i) > min(min_postforward))
            postforward_buf(i) = min(min_postforward);
        end
        if (prebackward_buf(i) > min(min_prebackward))
            prebackward_buf(i) = min(min_prebackward);
        end
        if (postbackward_buf(i) > min(min_postbackward))
            postbackward_buf(i) = min(min_postbackward);
        end
    end
end

