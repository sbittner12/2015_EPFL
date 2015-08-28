function h= plotCoherenceGroup(cohtype, subjects, exptype, channels, ...
                               musclestrs, trialstr, extraction_type)
%PLOTCOHERENCE Summary of this function goes here
%   plotCoherence(cohtype, subjects, exptype, channels, musclestr, trialstr, extraction_type)
%   Detailed explanation goes here
    %plottype = 'raw'; % ('raw', or 'sig', or raw_overlapping_sig', 'raw_multiply_sig')
nmuscles = length(musclestrs);
nchannelgroups = length(channels);
fontsize = 36;
save = 0;
format = 2;


switch format
    case 0
        figure;
        for j=1:nchannelgroups
            channel = channels{j};
            subplot(nmuscles+1, nchannelgroups+1, j+1);
            h1 = plot(0:1,0:1);
            ax = gca;
            str = '$$\longleftrightarrow$$';
            text(-.01, .5, sprintf('%s %s %s', channel), 'FontSize', fontsize, ...
                 'Interpreter', 'latex');
            set([h1,ax],'visible','off');
        end


        for i=1:nmuscles
            musclestr = musclestrs{i};
            fprintf('Plotting muscle: %s\n', musclestr);
            subplot(nmuscles+1, nchannelgroups+1, (i)*(nchannelgroups+1)+1)
            h1 = plot(0:1,0:1);
            ax = gca;
            str = '$$\longleftrightarrow$$';
            text(-.5, .5, sprintf('%s %s %s', musclestr), 'FontSize', fontsize, ...
                 'Interpreter', 'latex');
            set([h1,ax],'visible','off') 
            for j=1:nchannelgroups
                channel = channels{j};
                subplot(nmuscles+1, nchannelgroups+1, (i)*(nchannelgroups+1)+j+1), plotCoherence(cohtype, subjects, exptype, ...
                                           channel, musclestr, trialstr, format, extraction_type);
            end
        end
    case 2
        figure;
        if (length(channels) > 1)
            fprintf('Error: In format 2, must have only one channel.\n');
            return;
        end
        for i=1:nmuscles
            musclestr = musclestrs{i};
            fprintf('Plotting muscle: %s\n', musclestr);
            subplot(ceil(nmuscles / 5), 5, i) , plotCoherence(cohtype, subjects, exptype, ...
                                           channels{1}, musclestr, trialstr, format, extraction_type);
        end
        legend('Preprocessed EMG Signal');
end
