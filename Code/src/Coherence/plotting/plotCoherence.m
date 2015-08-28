function h= plotCoherence(cohtype, subjects, exptype, channels, musclestr, trialstr, format, extraction_type)
%PLOTCOHERENCE Summary of this function goes here
%   plotCoherence(cohtype, subjects, exptype, channels, musclestr, trialstr, extraction_type)
%   Detailed explanation goes here
    %plottype = 'raw'; % ('raw', or 'sig', or raw_overlapping_sig', 'raw_multiply_sig')
    max_ind_show = 26;
    freq_to_plot = [10, 20, 30, 40, 50];
    Fs = 1000;
    nfft = 512;   
    freq_res = Fs / nfft;
    fontsize = 24;
    fontsize2 = 16;
    fontsize3 = 10;
    nsubjects = length(subjects);
    resultsDir = ['E:',filesep,'Sean', filesep, 'Results', filesep];
    extractionstr = getExtractionstr(extraction_type);
    SUBJECTS = [3,4,5,6,7,8,9];
    plottype = 'raw';
    usebands = 0;
    bandlims = [4,8,14,20,30,50];
    bandlabels = {{'theta'; '(4-8Hz)'}, {'alpha'; '(8-14Hz)'}, {'beta'; '(14-20Hz)'}, ...
                  {'low gamma'; '(20-30Hz)'}, {'high gamma'; '(30-50Hz)'}};
    
    for i=1:nsubjects
        subject = subjects(i);
        cohDir = [resultsDir, 'Coherence', filesep, extractionstr, filesep, ...
                  sprintf('Subject%d', subject), filesep, 'RawCoherence', ...
                  filesep, musclestr, filesep];
        preEMGDir = [resultsDir, 'PreprocessedEMG', filesep, ...
                     sprintf('Subject%d', subject), filesep];
        sigDir = [resultsDir, 'Coherence', filesep, extractionstr, filesep, ...
                  sprintf('Subject%d', subject), 'SignificanceBootstrapping', filesep];
        if (exist(cohDir, 'dir') ~= 7)
            mkdir(cohDir);
        end
          
        mapfname = [cohDir, sprintf('Coherences_%s_%s_%s.mat', exptype, channels, trialstr)];
        sigfname = [sigDir, sprintf('SignificanceValues_%s_%s_%s_%s.mat', exptype, musclestr, channels, trialstr)];
        preEMGfname = [preEMGDir, sprintf('EMGPreproc_%s_%s_%s.mat', exptype, musclestr, trialstr)];

        load(mapfname);
        load(preEMGfname);
        [map, sig, cohname] = getCohmap(cohtype, mapfname, sigfname);
        Sig = (map - sig) > 0;      

        if (i == 1)
            freq_samples = size(map,1);
            time_samples = size(map,2);
            M_all = zeros(freq_samples, time_samples);
            Sig_all = zeros(freq_samples, time_samples);
            PreprocEMG = Preproc_percentile;
        else
            PreprocEMG = PreprocEMG + Preproc_percentile;
        end
        switch plottype
            case 'sig'
                Sig_all = Sig_all + Sig;
            case 'raw'
                M_all = M_all + map;
            case 'raw_overlapping_sig'
                map((map - sig) < 0) = 0;
                M_all = M_all + map;
            case 'raw_multiply_sig'
                Sig_all = Sig_all + Sig;
                M_all = M_all + map;
        end
        
        
%         emgdata = extractEMGSignal(subject, exptype, musclestr, trial);
%         emg = extractEMGSignalFast(emgdata, Event, targets, 'reach duration');
    end
    
    Sig_all = Sig_all / nsubjects;
    M_all = M_all / nsubjects;     
    PreprocEMG = PreprocEMG / nsubjects;
    
    if (nsubjects > 1)
        resDir = [resultsDir, 'Coherence', filesep, extractionstr, filesep, ...
                  'AllSubjects', filesep, 'plots', filesep, cohtype, filesep, ...
                  plottype, filesep];
    else
        resDir = [resultsDir, 'Coherence', filesep, extractionstr, filesep, ...
                  sprintf('Subject%d', subject), filesep, 'plots', filesep, ...
                  cohtype, filesep, plottype, filesep];
    end
    
    if (exist(resDir, 'dir') ~= 7)
            mkdir(resDir);
    end
    resfname = [resDir, sprintf('%s_%s_%s_%s_%s', exptype, musclestr, channels, trialstr, plottype)];
    
    sampling_windows = size(M_all,2);
    
    [cohmin, cohmax] = getCoherenceExtrema(subject, cohtype, extraction_type, trialstr);
    
    if (format == 1)
        h = figure;
        subplot(6,1,1);
        h1 = plot(0:1,0:1);
        ax = gca;
        str = '$$\longleftrightarrow$$';
        text(-.01, .5, sprintf('%s %s %s', musclestr, str, channels), 'FontSize', 42, ...
             'Interpreter', 'latex');
        set([h1,ax],'visible','off') 
        subplot(6,1,2:6), 
    end
    
    max_time_ind = nnz(~isnan(M_all(1,:)));
    
    MxSig_all = M_all.*Sig_all;
    %% Band plotting
    if (usebands)
        Sig_all = averageCoherenceBands(Sig_all, bandlims, freq_res);
        M_all = averageCoherenceBands(M_all, bandlims, freq_res);
        MxSig_all = averageCoherenceBands(MxSig_all, bandlims, freq_res);
        switch plottype
            case 'sig'
                X = Sig_all(:,1:max_time_ind);
            case 'raw'
                X = M_all(:,1:max_time_ind);
            case 'raw_overlapping_sig'
                %fprintf('muscle: %s %f\n', musclestr, max(max(M_all)));
                X = M_all(:,1:max_time_ind);
            case 'raw_multiply_sig'
                X = MxSig_all(:,1:max_time_ind);
        end
        maxX = max(max(X));
        buf = .05*maxX;
        
        [samples, windows] = size(PreprocEMG);
        nsamples_p = samples*windows;
        A = .02/(max(max(PreprocEMG)) - min(min(PreprocEMG)));
        emg = A*reshape(PreprocEMG, nsamples_p,1);
        meanoffset1 = mean(emg);
        yoffset1=maxX/2;
        
        figure;
        nbands = (length(bandlims)-1);
        for i=1:nbands;
            bandind = nbands-(i-1);
            subplot(nbands,1,i), plot(linspace(.5,11.5, nsamples_p), emg-meanoffset1+yoffset1, 'r');
            hold on;
            plot(1:11, X(bandind,:), 'b-*');
            ylabel(bandlabels{bandind}, 'FontSize', fontsize2);
            axis([1, max_time_ind, 0, maxX+buf]);
            ax = gca;
            set(ax,'XTick',1:2:11);
            if (i==length(bandlims)-1)
                set(ax,'XTickLabel',{'', '20%', '40%', '60%', '80%', '100%'}, 'FontSize', fontsize2);
                legend('EMG', 'Coherence');
            else
                set(ax,'XTickLabel',{''}, 'FontSize', fontsize2);
                hold off;
            end
        end
        return
     %% Not bands  
        
        
        
        
    else
        switch plottype
            case 'sig'
                imagesc(Sig_all(1:max_ind_show,1:max_time_ind), [0,1]);
            case 'raw'

                imagesc(M_all(1:max_ind_show,1:max_time_ind), [cohmin, cohmax]);
            case 'raw_overlapping_sig'
                %fprintf('muscle: %s %f\n', musclestr, max(max(M_all)));
                imagesc(M_all(1:max_ind_show,1:max_time_ind), [cohmin, cohmax]);
            case 'raw_multiply_sig'
                imagesc(MxSig_all(1:max_ind_show,1:max_time_ind), [cohmin, cohmax]);
        end

        if (format == 2)
            title(musclestr, 'FontSize', 40);
        end

        set(gca,'YDir','normal');
        %% Setup axes
        % Y axis
        ax = gca;
        freq_res = Fs/nfft;
        if (usebands)
            yaxis_pos = (freq_to_plot / bandwidth);
        else
            yaxis_pos = (freq_to_plot / freq_res);
        end
        set(ax,'YTick',yaxis_pos);
        % make tick labels
        nlabels = length(freq_to_plot);
        yaxis_labels = cell(1,nlabels); 
        for i=1:nlabels
            yaxis_labels{i} = sprintf('%d Hz', freq_to_plot(i));
        end
        set(ax,'YTickLabel',yaxis_labels, 'FontSize', fontsize);

        hold on;

        % Plot the start and uniform length buffer sizes
        if (strcmp(extraction_type, 'uniform length'))
            duration = premove_buf + postmove_buf_actual;
            bordershiftl = .5; % used to place trigger at beginning of plot
            bordershiftr = .5; % used to place end2 at end
            if (max_ind_show == 26)
                vshift = 2;
            else
                vshift = 10;
            end       
            begin_pos = 1-bordershiftl;
            end_pos = sampling_windows+bordershiftr;
            start_pos = sampling_windows*(premove_buf/duration);
            % X axis
            set(ax, 'XTick', [begin_pos, end_pos]);
            xaxis_labels{1} = sprintf('-%.2fs', premove_buf);
            xaxis_labels{2} = sprintf('+%.2fs', postmove_buf_actual);
            set(ax, 'XTickLabel',xaxis_labels);

            linewidth = 3;
            green_color = [34,139,34]/255;
            plot([start_pos, start_pos], [0, max_ind_show+1], 'Color', green_color, 'LineWidth', linewidth);
            %text(.5, max_ind_show+vshift, 'move onset', 'FontSize', fontsize);
            xlabel('time (s)');
            ylabel('Frequency (Hz)');
        elseif (strcmp(extraction_type, 'percentiles'))
            set(ax,'XTick',1:2:11);
            if (format == 1)
                set(ax,'XTickLabel',{'0%', '20%', '40%', '60%', '80%', '100%'}, 'FontSize', fontsize);
            else
                set(ax,'XTickLabel',{'0%', '20%', '40%', '60%', '80%', '100%'}, 'FontSize', fontsize-10);
            end
        end

        switch extraction_type
            case 'percentiles'
                [samples, windows] = size(PreprocEMG);
                nsamples_p = samples*windows;
                A = 5/(max(max(PreprocEMG)) - min(min(PreprocEMG)));
                emg = A*reshape(PreprocEMG, nsamples_p,1);
                meanoffset1 = mean(emg);
                yoffset1=5;
                plot(linspace(.5,11.5, nsamples_p), emg-meanoffset1+yoffset1, 'r');
            case 'reach duration'
                samples = sum(~isnan(map(1,:)));
                A = 5/(max(max(Preproc_reach_duration)) - min(min(Preproc_reach_duration)));
                emg = A*Preproc_reach_duration;
                nsamples_p = length(emg);
                meanoffset1 = mean(emg);
                yoffset1=5;
                plot(linspace(.5,samples+.5, nsamples_p), emg-meanoffset1+yoffset1, 'r');
        end
    end


    
    if (format == 1)
        saveas(h, resfname, 'png');
    end

