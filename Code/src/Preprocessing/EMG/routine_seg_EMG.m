function routine_seg_EMG

clear all;

Sub = [{'3'},{'4'},{'5'},{'6'},{'7'},{'8'},{'9'},{'10'}];
cond_EMG = [{'OAF_'},{'IAF_'},{'EE_healthy_'},{'EE_linear_'},{'EE_MJ_'}];
name = [{'Dominque'},{'Thomas'},{'Luca'},{'Matteo'},{'Miroslav'},{'Nicolas'},{'Beryl'},{'Stefano'}];
cond_Event = [{'OAF_'},{'IAF_'},{'HM_'},{'LT_'},{'MJ_'}];  

% correct_block = ['2','3','6','7','9']; %to change everytime we change of condition -  only put 5

%% Routine
kk = 1;
contin=1;
for s = 1:8
    for cc = 1:5
        if cc == 1
            n_block=12;
        else
            n_block=10;
        end
        for ii = 1:n_block; %block
            cd(['E:\Aurelie\Data\Event\Subject_' Sub{1,s} '\'])
            file_ev = [ cond_Event{1,cc} num2str(ii) '_event.mat'];
            if exist(file_ev, 'file') == 2
                load([ cond_Event{1,cc} num2str(ii) '_event.mat']);%change condition row every condition
                EMG = getfield(Event,'EMG');
                Start = getfield(EMG, 'Start');
                End = getfield(EMG,'End');
                Mov = zeros(16,1);
                Mov(1:2:end,1) = Start(1:2:end,1);
                Mov(2:2:end,1) = End(2:2:end,1);
                
                cd(['E:\Aurelie\Data\EMG\Raw_Data\Subject' Sub{1,s} '\' name{1,s} '\C3d\'])
                file_em = [cond_EMG{1,cc} num2str(ii) '_p.mat'];
                if exist(file_em, 'file') == 2
                    load(['E:\Aurelie\Data\EMG\Raw_Data\Subject' Sub{1,s} '\' name{1,s} '\C3d\' cond_EMG{1,cc} num2str(ii) '_p.mat'])
                    %             EMG_data = getfield(EMG,'data');
                    EMG_data = EMG_proc;
                    
                    for nn = 1:size(Mov,1)
                        if Mov(nn,:) > size(EMG_data,1)
                            Rejected_files{kk} = ['Sub' Sub{1,s} '_' cond_Event{1,cc} num2str(ii) ]
                            kk=kk+1;
                            contin=0;
                            break
                        end
                    end
                    if contin == 1
                        M1 = EMG_data(Mov(1,1):Mov(2,1),:); % all muscle
                        M2 = EMG_data(Mov(3,1):Mov(4,1),:);
                        M3 = EMG_data(Mov(5,1):Mov(6,1),:);
                        M4 = EMG_data(Mov(7,1):Mov(8,1),:);
                        M5 = EMG_data(Mov(9,1):Mov(10,1),:);
                        M6 = EMG_data(Mov(11,1):Mov(12,1),:);
                        M7 = EMG_data(Mov(13,1):Mov(14,1),:);
                        M8 = EMG_data(Mov(15,1):Mov(16,1),:);
                        
                        filename = ['E:\Aurelie\Data\Segmentation\ALL\Sub' Sub{1,s} '_' cond_Event{1,cc} '_' num2str(ii) ];
                        save(filename, 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8')
                    end
                else
                    Rejected_files{kk} = ['Sub' Sub{1,s} '_' cond_Event{1,cc} num2str(ii) ]
                    kk=kk+1;
                end
            else
                Rejected_files{kk} = ['Sub' Sub{1,s} '_' cond_Event{1,cc} num2str(ii) ]
                kk=kk+1;
            end
        end
        contin=1;
        clear M1 M2 M3 M4 M5 M6 M7 M8 Mov EMG EMG_data
    end
end
end


  