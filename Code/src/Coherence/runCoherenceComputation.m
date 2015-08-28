clear all
close all
eeggroups = {'parietal', 'center', 'frontal', 'left', 'right'};
muscles = {'PEC', 'INFRA', 'LAT', 'RHO', 'BICL', 'BICS' 'TRILAT', 'TRILONG', 'BRAC', 'PRO'};
nmuscles = length(muscles);
neeggroups = length(eeggroups);
subject = 3;
for i=1:nmuscles
    muscle = muscles{i};
    fprintf('Muscle: %s\n', muscle);
    for j=1:neeggroups
        eeggroup = eeggroups{j};
        fprintf('EEG %s\n', eeggroup);
        computeCoherences(subject, 'OAF', muscle, eeggroup, 1:12, 1);
    end
end