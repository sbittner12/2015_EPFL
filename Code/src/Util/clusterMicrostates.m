function [microstates_cl, ref_maps, groups]=clusterMicrostates(microstates, subjects, plot_flag)
%%CLUSTERMICROSTATES - finds corresponding microstates for each subject
%   clusterMicrostates(microstates, subjects, plot)
%
%   An ordering of the K microstates for each subject is found so that
%   corresponding microstates across subjects belong to the same component.
%   This is done by using k-means clustering on all microstates with K
%   clusters.  Since k-means does not guarantee that all clusters will have
%   the same number of elements, we rearrange the membership of the
%   microstates, so that each subject has one microstate belonging to one
%   cluster.  Cluster reference maps are determined by selecting the  
%   microstate of each group that has the maximum dot product with the rest
%   of the cluster.
%
%   PARAMETERS:
%     microstates - cell (1 x nsubjects) of matrices (nchans x K) : 
%                   the microstates for each subject
%     subjects     - vector (nsubjects) : the indices of each subject
%     plot        - scalar (1 or 0) : plot microstates before and after
%                   clustering
%   RETURNS:
%     microstates_cl - cell(1 x nsubjects) of matrices (nchans x K):
%                   the ordered microstates of each subject post clustering
%     ref_maps    - vector (K) : index of cluster reference maps
%                    Ex: For 2nd subject, 4th microstate, ind = 2*K + 4
%
%
%   @ 2015 Elvira Pirondini elvira.pirondini@epfl.ch (original author)
%   
%   @ 2015 Sean Bittner sbittner@andrew.cmu.edu      (adapted for project)
%

nsubjects = length(subjects);
K = size(microstates{1},2);
chanlocs = loadChanlocs('all');
nchans = length(chanlocs);

if (plot_flag)
    plotMicrostates(microstates, subjects, chanlocs);
    suptitle('Microstates before clustering');
end

% Setup struct array of microstate information for clustering algorithm
mins = zeros(1,nsubjects);
maxs = zeros(1,nsubjects);
for i=1:nsubjects
    mins(i) = min(min(microstates{i}));
    maxs(i) = max(max(microstates{i}));
    for k=1:K
        reaching_microstates(K*(i-1)+k).microstate = microstates{i}(:,k);
        reaching_microstates(K*(i-1)+k).type = '';
        reaching_microstates(K*(i-1)+k).subject = i;
        reaching_microstates(K*(i-1)+k).chanlocs = chanlocs;
        reaching_microstates(K*(i-1)+k).nchans = nchans;
    end
end
% Used for scaling microstate maps from topoplot
mini = min(mins);
maxi = max(maxs);
    
% Fill x (ntotalmics x n) which holds rows of interpolated microstate maps
ntotalmics=length(reaching_microstates);
subjects2 = zeros(ntotalmics,1);
for ii=1:ntotalmics
    [~,map,~,~,~]=topoplot(reaching_microstates(ii).microstate,reaching_microstates(ii).chanlocs,...
        'conv','on','electrodes','off','style','map','noplot','on','maplimits',[mini,maxi]);
    if (ii==1)
        n = size(map(~isnan(map)), 1);
        x = zeros(ntotalmics, n);
    end
    x(ii,:)=map(~isnan(map));
    subjects2(ii)=reaching_microstates(ii).subject;
end

% Clustering
[clusters,ref_maps,groups]=clustering1(x,K,subjects2);

[microstates_cl, corr_cl] = getClusteredMicrostates(microstates, K, groups);


if (plot_flag)
    plotMicrostates(microstates_cl, subjects, chanlocs, corr_cl);
    suptitle('Microstates after clustering');
end

