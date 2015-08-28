function [h] = plotEvent(Data)
%PLOTEVENT Summary of this function goes here
%   Detailed explanation goes here
h = figure;
plot(1:8, Data.Trigger, '*', 1:8, Data.Start(1:2:15), '*', ...
     1:8,  Data.Start(2:2:16), '*', 1:8, Data.End(1:2:15), '*', ...
     1:8, Data.End(2:2:16), '*');
legend('trigger', 'start1', 'start2', 'end1', 'end2');

end

