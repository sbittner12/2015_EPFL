%% Function for clustering
subjects = [3,4,5];
nsubjects = length(subjects);
nmicrostates = 6;

% maps: are the topographical maps of each subject
EEG = loadEEG(3,'OAF','all');
chanlocs = EEG(1).chanlocs;
nchans = length(chanlocs);

microstates = loadMicrostates(subjects, nmicrostates);
plotMicrostates(microstates, subjects, chanlocs);

mins = zeros(1,nsubjects);
maxs = zeros(1,nsubjects);
for i=1:nsubjects
    mins(i) = min(min(microstates{i}));
    maxs(i) = max(max(microstates{i}));
    for j=1:nmicrostates
        reaching_microstates(nmicrostates*(i-1)+j).microstate = microstates{i}(:,j);
        reaching_microstates(nmicrostates*(i-1)+j).type = '';
        reaching_microstates(nmicrostates*(i-1)+j).subject = i;
        reaching_microstates(nmicrostates*(i-1)+j).chanlocs = chanlocs;
        reaching_microstates(nmicrostates*(i-1)+j).nchans = nchans;
    end
end
mini = min(mins);
maxi = max(maxs);
    

%reaching_microstates=struct([reaching_microstates struct('microstate',maps(:,nn),'type',type,'subject',ii,'chanlocs',chanlocs,'nchans',nchans)]);

n_mics=length(reaching_microstates);
subjects2 = zeros(n_mics,1);
for ii=1:n_mics
    [~,map,~,~,~]=topoplot(reaching_microstates(ii).microstate,reaching_microstates(ii).chanlocs,...
        'conv','on','electrodes','off','style','map','noplot','on','maplimits',[mini,maxi]);
    if (ii==1)
        n = size(map(~isnan(map)), 1);
        x = zeros(n_mics, n);
    end
    x(ii,:)=map(~isnan(map));
    subjects2(ii)=reaching_microstates(ii).subject;
end

% Clustering
[clusters,ref_maps,groups]=clustering1(x,nmicrostates,subjects2);

microstates_clustered = getClusteredMicrostates(microstates, nsubjects, nmicrostates, groups);

plotMicrostates(microstates_clustered, subjects, chanlocs);

