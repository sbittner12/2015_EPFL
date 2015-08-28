function [] = plotMicrostates(microstates, subjects, chanlocs, varargin)
%PLOTMICROSTATES Summary of this function goes here
%   Detailed explanation goes here
    figure;
    corr = 0;
    if (length(varargin) > 0)
        corr = 1;
        corr_cl = varargin{1};
    end
    
    nsubjects = length(subjects);
    nmicrostates = size(microstates{1},2);
    for i=1:nsubjects
        subject = subjects(i);
        for j=1:nmicrostates
            ind = nmicrostates*(i-1)+j;
            subplot(nsubjects, nmicrostates, ind);
            topoplot(microstates{i}(:,j), chanlocs);
            if (corr)
                if (j==1)
                    ylabel(sprintf('Subject %d', subject));
                    title(sprintf('%d: %f', j, corr_cl{i}(j)));
                else
                    title(sprintf('%d: %f', j, corr_cl{i}(j)));
                end
            else
                if (j==1)
                    ylabel(sprintf('Subject %d', subject));
                    title(sprintf('%d', j));
                else
                    title(sprintf('%d', j));
                end
            end
        end
    end

end

