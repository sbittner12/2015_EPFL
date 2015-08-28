function [microstates_discrete] = discretizeMicrostates(microstate_activations, nsamples)
%DISCRETIZEMICROSTATES Summary of this function goes here
%   Detailed explanation goes here
    microstates_discrete = zeros(size(microstate_activations));
    for i=1:nsamples
        [x, ind] = max(microstate_activations(i,:));
        microstates_discrete(i,ind) = 1;
    end
end

