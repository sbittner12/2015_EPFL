function [gPDC] = computegPDC(C, A_w, nfreqs)
    % Calculate generalized PDC
    nv = size(C,1);
    nfftsamples = nfreqs/2+1;
    gPDC = zeros(nv,nv,nfftsamples);
    for j=1:nv
        for freq_ind=2:nfftsamples
            gPDC_num = zeros(nv,1);
            gPDC_denom = 0;
            for i=1:nv
                gPDC_num(i) = abs(A_w{freq_ind}(i,j))^2 / C(i,i);
                gPDC_denom = gPDC_denom + gPDC_num(i);
            end
            gPDC(:,j, freq_ind) = gPDC_num / gPDC_denom;
        end
    end

end