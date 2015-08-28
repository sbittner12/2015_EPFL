function [target_begin_ind, target_end_ind] = parseTargets(target, Event)
%PARSETARGETS Parses the target input parameter
%   parseTargets(target, Event) returns the target indices associated with 
%   the target value.  If the target value is a target number, it will use 
%   the corresponding Event structure to determine the index based on the
%   target order.  If the target parameter is 'all'.  It will return a 1
%   for the begining index and an 8 for the final index.
%
%  PARAMETERS
%    target - integer specifying a single target or the string 'all' 
%             specifying that all targets are desired.
%    Event - event structure containing the target order.
%
%  RETURNS
%    target_begin_ind - index of first target
%    target_end_ind - index of final target (same as beginning for single
%                                            target)
%
%  @ 2015 Sean Bittner   sbittner@andrew.cmu.edu
%

    if (strcmp(target, 'all'))
        target_begin_ind = 1;
        target_end_ind = 8;
    else
        target_begin_ind = find((Event.Target_order == target), 1);
        target_end_ind = target_begin_ind;
    end
end

