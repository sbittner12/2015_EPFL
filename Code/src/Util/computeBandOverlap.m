function [overlap] = computeBandOverlap(s1,e1,s2,e2)
%COMPUTEBANDOVERLAP Summary of this function goes here
%   Detailed explanation goes here
if (e1 <= e2)
    if (s1 >= s2)
        overlap = e1-s1;
    else
        overlap = e1-s2;
    end
else
    if (s1 >= s2)
        overlap = e2-s1;
    else
        overlap = e2-s2;
    end

end

