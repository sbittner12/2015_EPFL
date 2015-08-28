function [clusters,ref_maps,groups]=clustering1(x, K, subjects)
%%CLUSTERING1 - clusters data, each cluster has one for each subject
%   clustering1(x, n, subjects)
%
%   Each row of x is clustered sorted using k-means clustering with n
%   clusters.  After this preliminary clustering assignment, we re-cluster
%   the data so that each cluster has one data element from each subject.
%
%   PARAMETERS:
%     x        - matrix (nsubjects*K x n) : data elements
%     K        - scalar > 0 : number of clusters for k-means
%     subjects - vector (nsubjects) : the indices of each subject

%   RETURNS:
%     clusters - cell(K x 1) of vectors : the indices of the data elements
%                that belong to each cluster from k-means clustering
%     ref_maps - vector (K) : the indices of the data elements that are
%                                most representative of each cluster.
%     groups   - cell (K x 1) of matrices (nsubjects x 2): for each
%                microstate, the index of each subject's member microstate
%                is in the first column, and the second column is the
%                microstate's correlation with the rest of the group
%
%   @ 2015 Elvira Pirondini elvira.pirondini@epfl.ch (original author)
%   
%   @ 2015 Sean Bittner sbittner@andrew.cmu.edu      (adapted for project)
%

% Global variables
nsubjects=length(unique(subjects));
% for ii=1:size(x,1)
%     x(ii,:)=x(ii,:)-x(ii,1);
% end
clust=kmeans(abs(x),K,'distance','correlation', 'replicates', 50);
for ii=1:K
    clusters{ii}=find(clust==ii);
    len=length(clusters{ii});
    corr_coef_max=0;
    for jj=1:len
        corr_coef=0;
        for kk=1:len
            corr_coef=corr_coef+abs(corr2(x(clusters{ii}(jj),:),...
                x(clusters{ii}(kk),:)));
        end
        if(corr_coef>corr_coef_max)
            corr_coef_max=corr_coef;
            ref_map=clusters{ii}(jj);
        end
    end
    ref_maps(ii)=ref_map;
end
% ref_maps=sort(ref_maps);
groups=cell(1,K);
map1=0;
% For each subject
for ii=1:nsubjects
%     subject = subjects(ii);
    nmics=sum(subjects==ii);
    co=zeros(K,nmics);
    % For each microstate
    for jj=1:nmics
        map1=map1+1;
        for kk=1:K
            map2=ref_maps(kk);
            if(ii~=subjects(map2))
                co(kk,jj)=abs(corr2(x(map1,:),x(map2,:)));
            elseif(map1==map2)
                co(kk,jj)=1;
            else
                co(kk,jj)=0;
            end
        end
    end
    for jj=1:nmics
        [maxs,inds]=max(co);
        [max_co,col]=max(maxs);
        row=inds(col);
        groups{row}=[groups{row};[map1-nmics+col,max_co]];
        co(row,:)=0;
        co(:,col)=0;
    end
end