function [microstates_cl, corr_cl] = getClusteredMicrostates(microstates, K, groups)
%%GETCLUSTEREDMICROSTATES - decode group indices, and extract from original
%   getClusteredMicrostates(x, n, subjects)
%  
%   Creates a cell of subject microstates according to the group
%   specification from clustering.
%  
%   PARAMETERS:
%     microstates - cell (1 x nsubjects) of matrices (nchans x K) : 
%                   the microstates for each subject
%     K           - scalar > 0 : number of clusters for k-means
%     groups      - cell (K x 1) of matrices (nsubjects x 2): for each
%                   microstate, the index of each subject's member micro-
%                   state is in the first column, and the second column is 
%                   the microstate's correlation with the rest of the group
%   RETURNS:
%     microstates_cl - cell(1 x nsubjects) of matrices (nchans x K):
%                   the ordered microstates of each subject post clustering
%   
%   @ 2015 Sean Bittner sbittner@andrew.cmu.edu
%
    nsubjects = length(microstates);
    microstates_cl = cell(nsubjects,1);
    for i=1:nsubjects
        for j=1:K
            ind = groups{j}(i,1);
            ii = floor((ind-1) / K) + 1;
            jj = mod((ind-1), K)+1;
            microstates_cl{i}(:,j) = microstates{ii}(:,jj);
            corr_cl{i}(j) = groups{j}(i,2);
        end
    end
end

