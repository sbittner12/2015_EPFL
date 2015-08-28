function [begin_ind, end_ind] = getReachIndsReachDuration(Data, Event, targets, srate, srate_new)
%GETREACHINDSUNIFORMLENGTH Summary of this function goes here
%   Detailed explanation goes here
    [begin_target_ind, end_target_ind] = parseTargets(targets, Event);
    trigger = parseTriggers(Data.Trigger, Data.End(end));
    begin_ind = round(trigger(begin_target_ind)*srate_new/srate);
    end_ind = round(Data.End(2*end_target_ind)*srate_new/srate);
end

