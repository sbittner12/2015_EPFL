function [premove_buf, postmove_buf] = getMovementBuffers(subjects, varargin)
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

premove_buf = 100*ones(nexptypes,1);
postmove_buf = 100*ones(nexptypes,1);

for k = 1:nsubjects
    subject = subjects(k);
    subjectEventDir = [dataDir, 'Event', filesep, sprintf('Subject_%d', subject), filesep];
    
    for i=1:nexptypes
        exptype = exptypes{i};
        %fprintf('******** %s ********\n', exptype);
        max_premovement = [];
        min_premovement = [];
        max_postmovement = [];
        min_postmovement = [];
        ind = 1;
        for j=1:ntrials
            eventfname = [subjectEventDir, sprintf('%s_%d_event.mat',exptype, j)];
            if (exist(eventfname, 'file') == 2)
                load(eventfname);
                triggers = parseTriggers(Event.EEG.Trigger, Event.EEG.End(end));
                % Find min and max premovement periods
                premovement = (Event.EEG.Start(1:2:end-1) - triggers) / EEG_Fs;
                max_premovement(ind) = max(premovement);
                min_premovement(ind) = min(premovement);
                % Find min and max postmovement periods
                postmovement = (Event.EEG.End(2:2:end) - Event.EEG.Start(1:2:end-1)) / EEG_Fs;
                max_postmovement(ind) = max(postmovement);
                min_postmovement(ind) = min(postmovement);
                ind = ind + 1;
            else
                %display('file dne');
            end
        end
        if (premove_buf(i) > min(min_premovement))
            premove_buf(i) = min(min_premovement);
        end
        if (postmove_buf(i) > min(min_postmovement))
            postmove_buf(i) = min(min_postmovement);
        end
%         fprintf('Max pre-movement period :%.2f\n', max(max_premovement));
%         fprintf('Min pre-movement period :%.2f\n', min(min_premovement));
%         fprintf('Max post-movement period :%.2f\n', max(max_postmovement));
%         fprintf('Min post-movement period :%.2f\n', min(min_postmovement));
    end
end

