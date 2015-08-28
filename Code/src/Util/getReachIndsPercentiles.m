function [begin_ind, end_ind] = getReachIndsPerentiles(Data, Event, targets, srate, srate_new, windowsize)
%GETREACHINDSUNIFORMLENGTH Summary of this function goes here
%   Detailed explanation goes here
    [begin_target_ind, end_target_ind] = parseTargets(targets, Event);
    prestart_buf = (windowsize / 2) - 1;
    postend_buf = windowsize / 2;
    trigger = parseTriggers(Data.Trigger, Data.End(end));
    begin_ind = round(trigger(begin_target_ind)*srate_new/srate) - prestart_buf;
    end_ind = round(Data.End(2*end_target_ind)*srate_new/srate) + postend_buf;
end

