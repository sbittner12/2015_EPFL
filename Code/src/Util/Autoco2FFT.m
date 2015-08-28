function A_w = Autoco2FFT(acfunc, nfreqs)
%AUTOCO2FFT Summary of this function goes here
%   Detailed explanation goes here
p = length(acfunc) - 1;
nv = size(acfunc{1},1);
for i = 1:nv
    for j = 1:nv 
        vr = zeros(nfreqs,1); %fillchar1(vr,0);
        if i == j 
            vr(1) =1;
        end
    
        for t = 2 : p+1
          vr(t) = - acfunc{t}(i,j);
        end
        vi2 = fft(vr); %% questo è da controllare perchè non sono sicura!
       
        for t = 2:nfreqs 
            A_w{t}(i,j) = vi2(t);
        end
    end
end 

A_w{1} = zeros(nv);

end

