function [avgsig] = crossTrialAverage(x, trial_length, ntrials)
%%CROSSTRIALAVERAGE - averages signals across uniform length trials
%   crossTrialAverage(x, trial_length, ntrials)
%
%   PARAMETERS:
%     x            - vector (ntrials*trial_length) : signal
%     trial_length - scalar : uniform length of each trial
%     ntrials      - scalar : total number of trials in x
%
%   RETURNS:
%     avgsig       - vector (trial_length) : trial-averaged signal
%
%   @ 2015 Sean Bittner sbittner@andrew.cmu.edu   
%
assert(length(x) == ntrials*trial_length);
avgsig = mean(reshape(x, trial_length, 12)');

end

