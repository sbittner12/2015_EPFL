function X = getReachInds(Event, datatype, targets, srate, srate_new)
%GETREACHINDS Summary of this function goes here
%   Detailed explanation goes here 
    if (strcmp(datatype, 'EEG'))
        Data = Event.EEG;
    else
        Data = Event.EMG;
    end
    trigger = parseTriggers(Data.Trigger, Data.End(end));
    [target_begin_ind, target_end_ind] = parseTargets(targets, Event);
    trigger_ind = round(trigger(target_begin_ind)*srate_new/srate);
    start1_ind = round(Data.Start((2*target_begin_ind)-1)*srate_new/srate);
    end1_ind   = round(Data.End((2*target_begin_ind)-1)*srate_new/srate);
    start2_ind = round(Data.Start(2*target_end_ind)*srate_new/srate);
    end2_ind =   round(Data.End(2*target_end_ind)*srate_new/srate);

    X = [trigger_ind; start1_ind; end1_ind; start2_ind; end2_ind;];
end

