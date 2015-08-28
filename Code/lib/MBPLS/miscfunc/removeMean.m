function output = removeMean(data)

num_channels = size(data,1);

for i=1:num_channels
    data(i,:) = data(i,:)-mean(data(i,:));
end

output = data;