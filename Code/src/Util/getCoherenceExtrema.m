function [cohmin, cohmax] = getCoherenceExtrema(subjects, cohtype, extraction_type, trialstr)
%FINDCOHERENCEEXTREMA Summary of this function goes here
%   Detailed explanation goes here
exptype = 'OAF';
muscles = {'TRAPS', 'TRAPM', 'DANT', 'DMED', 'DPOS', 'PEC', 'INFRA', 'LAT', 'RHO', 'BICL', 'BICS' 'TRILAT', 'TRILONG', 'BRAC', 'PRO'};
extractionstr = getExtractionstr(extraction_type);
cohDir = ['E:',filesep,'Sean', filesep, 'Results', filesep, 'Coherence', ... 
          filesep, extractionstr filesep];
%eeggroups = {'frontal', 'parietal', 'left', 'right', 'center'};
eeggroups = {'parietal'};
neeggroups = length(eeggroups);
nsubjects = length(subjects);
nmuscles = length(muscles);
init_flag = 0;
for i=1:nsubjects
    subject = subjects(i);
    subjectDir = [cohDir, sprintf('Subject%d', subject), filesep];
    for j=1:neeggroups
        eeggroup = eeggroups{j};
        for k=1:nmuscles
            musclestr = muscles{k};
            mapfname = [subjectDir, 'RawCoherence', filesep, musclestr, filesep, ...
                        sprintf('Coherences_%s_%s_%s.mat', ...
                        exptype, eeggroup, trialstr)];
            sigfname = [subjectDir, 'SignificanceBootstrapping', filesep, ...
                        sprintf('SignificanceValues_%s_%s_%s_%s.mat', exptype, ...
                        musclestr, eeggroup, trialstr)];
            if (exist(sigfname, 'file') ~= 2)
                fprintf('DNE %s\n', sigfname);
                sigfname = '';
                if (exist(mapfname, 'file') ~= 2)
                    continue
                end
            else
                %fprintf('EXISTS %s\n', sigfname);
            end
            [map, sig, cohname] = getCohmap(cohtype, mapfname, sigfname);
            cohmin_ij = min(min(map));
            cohmax_ij = max(max(map));
            if (init_flag == 0)
                init_flag = 1;
                cohmin = cohmin_ij;
                cohmax = cohmax_ij;
            else
                if (cohmin_ij < cohmin)
                    cohmin = cohmin_ij;
                end
                if (cohmax_ij > cohmax)
                    cohmax = cohmax_ij;
                end
            end
        end
    end
end
    

