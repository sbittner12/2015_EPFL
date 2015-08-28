function [iCOH] = computeiCOH(C, A_w, nfreqs)
    nv = size(C,1);
    nfftsamples = nfreqs/2+1;
    iCOH = zeros(nv,nv,nfftsamples);
    for i=1:nv
        for j=1:nv
            if (i == j)
                iCOH(i,j,:) = zeros(nfftsamples,1);
            else
                for freq_ind=2:(nfftsamples)
                    iCOH_num = abs(A_w{freq_ind}(i,j))^2 / C(i,i);
                    iCOH_denom = iCOH_num + (abs(A_w{freq_ind}(j,j))^2) / C(j,j);
                    iCOH(i,j, freq_ind) = iCOH_num / iCOH_denom;
                end
            end
        end
    end
end

