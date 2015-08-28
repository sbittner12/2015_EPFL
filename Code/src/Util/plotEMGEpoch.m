function [h] = plotEMGEpoch(subject, exptype, trial, target)
%PLOTEMGEPOCH Summary of this function goes here
%   Detailed explanation goes here
    emgDir = ['E:' filesep, 'Sean', filesep, 'Data', filesep, 'EMG', filesep];
    emgEpochDir = [emgDir, sprintf('Subject%d', subject), filesep, ...
                   'SingleEpochMat', filesep];
    emgEpochfname = [emgEpochDir, sprintf('%s_%d_%d.mat', exptype, trial, target)];
    emgLabelsfname = [emgDir, 'EMG_labels.mat'];
    load(emgEpochfname);
    load(emgLabelsfname);    
    
    data = EMG.data;
    [Nemg, samples] = size(data);
    
    h = figure;
    for i=1:Nemg
        subplot(3,5,i), plot(linspace(0,1,samples), data(i,:));
        title(correct_labels{i});
    end
    
    suptitle(sprintf('Subject %d, Trial %d, Target %d', subject, trial, target));

end

