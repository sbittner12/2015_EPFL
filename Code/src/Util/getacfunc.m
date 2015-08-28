function acfunc = getacfunc(C, A)
%GETACFUNC Summary of this function goes here
%   Detailed explanation goes here
    nv = size(C,1);
    p = size(A,2)/2;
    acmat = [C, A];
    acfunc = cell(p*1, 1);
    for k =1:p+1
        for i = 1:nv
            for j = 1:nv
                acfunc{k}(i,j) = acmat(i,j+(k-1)*nv);
            end
        end
    end

end

