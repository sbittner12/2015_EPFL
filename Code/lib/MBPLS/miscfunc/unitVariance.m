function output = unitVariance(data)

[num_ch, N] = size(data);

for i = 1:num_ch
    sigma = std(data(i, :));
    output(i, :) = data(i, :)/sigma;
end