function [Qb_dotprods, Qb_colors] = rankQb(Qb, k)
%COLORQB Summary of this function goes here
%   Detailed explanation goes here
N = length(Qb);
nelectrodes = size(Qb{1},k);
bestfit_ind = 0;
bestfit_val = 0;
minfit_val = 0;
dotprods = zeros(N,1);

for i=1:N
    dotprod_sum = 0;
    for j=1:N
        if (i~=j)
            dotprod_sum = dotprod_sum + abs(Qb{i}(:,k))' * abs(Qb{j}(:,k));
        end
    end
    if (dotprod_sum > bestfit_val)
        bestfit_ind = i;
        bestfit_val = dotprod_sum;
    end
    if (dotprod_sum < minfit_val)
        minfit_val = dotprod_sum;
    end
    dotprods(i) = dotprod_sum;
end

for i=1:N
    if (i == bestfit_ind)
        Qb_colors{i} = [0,153,0]/255; % green
    else
        r = 255*(dotprods(i) - minfit_val) / (bestfit_val - minfit_val);
        g = 0;
        b = 255 - 255*(dotprods(i) - minfit_val) / (bestfit_val - minfit_val);
        Qb_colors{i} = [r,g,b]/255; 
    end
end

for i=1:N
    Qb_dotprods(i) = abs(Qb{i}(:,k))' * abs(Qb{bestfit_ind}(:,k));
end
Qb_dotprods = Qb_dotprods / max(Qb_dotprods);


