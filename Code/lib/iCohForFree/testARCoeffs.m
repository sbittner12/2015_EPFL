% Script for finding where AR coefficients fall within acmat variable of
% provided iCoh code
load pair1
p = 3;
nfft = 512;
[w, A, C, sbc, fpe, th] = arfit(pair1, p, p);
[bb, dim, acmat] = iCOH(pair1, p, nfft);
[rows, cols] = size(A);
dims = size(acmat, 1);
for i = 1:rows
    for j = 1:cols
        for ii = 1:dims
            for jj = (ii+1):dims
                th = abs(.02*A(i,j));
%                 fprintf('(%d,%d), -> (%d,%d)\n', i, j, ii, jj);
%                 fprintf('%f, -> %f\n\n', A(i,j), acmat(ii,jj));
                if (abs(A(i,j) - acmat(ii,jj)) < th)
                    fprintf('A(%d,%d) is close to acmat(%d,%d)\n', i, j, ii, jj);
                end
            end
        end
    end
end