function [trigger_out] = parseTriggers(triggers_in, last_end)
%PARSETRIGGERS Parses the Trigger property of an Event structure
%   parseTriggers(triggers_in, last_end) determines which of the elements
%   of the trigger vector are supposed to be kept.  There are only 8
%   targets and therefor epochs per trial, however, sometimes there are 9
%   or even 10 trigger indices for an Event structure.  I based this
%   trigger parsing on Aurelie's code from earlier this year.
%
%  PARAMETERS
%    trigger_in - integer vector of trigger indices
%    last_end - the last element of the end vector from the same Event.
%
%  RETURNS
%    trigger_out - a 1x8 vector of the parsed trigger indices
%
%  @ 2015 Sean Bittner   sbittner@andrew.cmu.edu
%

ntriggers = length(triggers_in);
if (ntriggers >= 10)
    trigger_out = triggers_in(2:end-1);
elseif (ntriggers == 9)
    if (triggers_in(end) > last_end)
        trigger_out = triggers_in(1:end-1);
    else
        trigger_out = triggers_in(2:end);
    end
else
    assert(ntriggers == 8);
    trigger_out = triggers_in;
end

end

