subjects = [4,5,7,8,9];
eeggroups = {'frontal', 'parietal', 'left', 'right', 'center'};
muscles = {'INFRA', 'LAT', 'RHO', 'BICL', 'BICS' 'TRILAT', 'TRILONG', 'BRAC', 'PRO'};
nsubjects = length(subjects);
nmuscles = length(muscles);
neeggroups = length(eeggroups);

for i=1:nmuscles
    musclestr = muscles{i};
    fprintf('************************* Muscle %s ***********************\n\n', musclestr);
    for k=1:nsubjects
        subject = subjects(k);
        fprintf('Subject %d\n', subject);
        for j=1:neeggroups
            eeggroup = eeggroups{j};
            fprintf('************************* EEG %s ***********************\n\n', eeggroup);
            computeSignificanceValues(subject, 'OAF', musclestr, eeggroup, 1:12,1:8);
        end
    end
end