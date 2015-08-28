function [x2perm] = Phase_randomization(x)

F = fft(x); % default fft operates by column: one regressor per column
FA = angle(F);
FM=abs(F);                
       
FA2 = FA([1 randperm(size(FA,1)-1)+1],:); % preserving DC un-permuted
F2 = FM.*exp(FA2*sqrt(-1));
x2perm = real(ifft(F2));

end