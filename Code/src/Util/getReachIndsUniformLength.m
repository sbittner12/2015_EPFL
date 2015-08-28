function [begin_ind, end_ind] = getReachIndsUniformLength(Data, Event, targets, srate, srate_new, premove_buf, postmove_buf)
%GETREACHINDSUNIFORMLENGTH Summary of this function goes here
%   Detailed explanation goes here
    [begin_target_ind, end_target_ind] = parseTargets(targets, Event);
    move_onset = Data.Start(2*begin_target_ind-1)*srate_new/srate;
    begin_ind = round(move_onset - (premove_buf*srate_new));
    end_ind = round(move_onset + (postmove_buf*srate_new));
end

