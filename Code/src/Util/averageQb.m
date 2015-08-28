function [Qbavg] = averageQb(Qb, comp)
%AVERAGEQB Summary of this function goes here
%   Detailed explanation goes here
nblocks = length(Qb);
[nmuscles, ncomps] = size(Qb{1});
assert(0 < comp  && comp <= ncomps);
Qbcomp = zeros(nblocks, nmuscles);

for i=1:nblocks
    Qbcomp(i,:) = abs(Qb{i}(:,comp))';
end

Qbavg = mean(Qbcomp);

end

