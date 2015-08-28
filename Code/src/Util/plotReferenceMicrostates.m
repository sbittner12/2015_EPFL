function [] = plotReferenceMicrostates(microstates, ref_maps, chanlocs)
%PLOTREFERENCEMICROSTATES Summary of this function goes here
%   Detailed explanation goes here
    nmicrostates = size(microstates{1},2);
    figure;
    for i=1:nmicrostates
        ind = ref_maps(i);
        ii = floor((ind-1) / nmicrostates) + 1;
        jj = mod((ind-1), nmicrostates) + 1;
        subplot(1,nmicrostates, i), topoplot(microstates{ii}(:,jj), chanlocs);
        title(sprintf('MS %d', i));
    end
    suptitle('Microstate reference patterns');
end