function [PDC] = computePDC(C, A_w, nfreqs)
    nv = size(C,1);
    nfftsamples = nfreqs/2+1;
    PDC = zeros(nv,nv,nfftsamples);
    for j=1:nv
        for freq_ind=2:(nfftsamples)
            PDC_num = zeros(nv,1);
            PDC_denom = 0;
            for i=1:nv
                PDC_num(i) = abs(A_w{freq_ind}(i,j))^2;
                PDC_denom = PDC_denom + PDC_num(i);
            end
            PDC(:,j, freq_ind) = PDC_num / PDC_denom;
        end
    end


end

