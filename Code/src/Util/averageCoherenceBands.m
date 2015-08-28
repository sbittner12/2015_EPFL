function [M_bands] = averageCoherenceBands(M, bandlims, freq_res)
%AVERAGECOHERENCEBANDS Summary of this function goes here
%   Detailed explanation goes here
    [freq_bins, time_bins] = size(M);
    nbands = length(bandlims) - 1;
    bandliminds = round(bandlims / freq_res);
    startinds = bandliminds(1:(end-1))+1;
    endinds = bandliminds(2:end);

    M_bands = zeros(nbands, time_bins);

    for i=1:nbands
        startind = startinds(i);
        endind = endinds(i);
        M_bands(i,:) = mean(M(startind:endind, :));
    end

end