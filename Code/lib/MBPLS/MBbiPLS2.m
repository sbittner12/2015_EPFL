function [Tb,Pb,Wb, Wb_reproj, Wt,Tt,Ub,Qb,Wu,Tu, ssx, ssy, Wridge] = MBbiPLS2(X,Y,nLV,tol)
% function [Tb,Pb,Wb,Wt,Tt,Ub,Qb,Wu,Tu,B,ssq,Rbo,Rbv,Ry,Lbo,Lbv,Lto] =
% MBbiPLS2(X,Y,nLV,tol)
% in  :
% X  (1 x number of X-blocks) cell array where each element is an X-data-block 
%   of dimension t (time) x  d (channel) 
% Y  (1 x number of Y-blocks) cell array where each element is an Y-data-block 
%   of dimension t (time) x d2(channel)
% nLV (1 x 1) number of latent variables extracted per X-block
% tol (1 x 1) tolerance for convergence (default 1e-12)
%
% out :
% Tb (objects x number of X-blocks.max(nLV)) block scores, [t1-block-1 t1-block-2 ...t2-block-1...]
% Pb (X-variables x max(nLV)) X-block loadings
% Wb (X-varaibles x max(nLV)) X-block weights
% Wt (number of X-blocks x max(nLV)) X-block super weights
% Tt (objects x max(nLV)) X-block super scores
% Ub (objects x number of Y-blocks.max(nLV)) Y-block scores
% Qb (Y-variables x max(nLV)) Y-block weights
% Wu (number of Y-blocks x max(nLV)) Y-block super weights
% Tu (objects x max(nLV)) Y-block super scores
% ssx
% ssy

if nargin == 0
   help MBbiPLS2
   return;
elseif nargin == 3
   tol = 1e-12;
end

maxiter = 10000;
nbX = length(X);


Xres = X;
Yres = Y;
Tb = cell(1, nbX); 
% Pb = cell(1, nbX); 
Pb = cell(1, nbX); 
Wb = cell(1, nbX); 
Wb_reproj = cell(1, nbX);
DimX = cell(1, nbX);
DimY = cell(1, nbX);
for b=1:nbX
    DimX{b} = size(X{b});
    DimY{b} = size(Y{b});
    Tb{b} = zeros(DimX{b}(1), nLV);
    Pb{b} = zeros(DimX{b}(2), nLV);
    Wb{b} = zeros(DimX{b}(2), nLV);
    Wb_reproj{b} = zeros(DimX{b}(2), nLV);
end
Wt = zeros(nbX,nLV);
Tt = zeros(DimX{1}(1), nLV);
Ub = cell(1, nbX); 
Qb = cell(1, nbX); 
for b = 1:nbX
    Ub{b} = zeros(DimY{b}(1), nLV);
    Qb{b} = zeros(DimY{b}(2), nLV);
end
Wu = zeros(nbX,nLV);
Tu = zeros(DimX{1}(1),nLV);
B = cell(1, nbX);
for b = 1:nbX
    B{b} = zeros(nLV,nLV);
end
Wridge = cell(1, nbX);

SSX = zeros(1, nbX);
SSY = zeros(1, nbX);
ssx = cell(1, nbX);
ssy = cell(1, nbX);
for aa = 1:nbX
    SSX(aa)=sum(sum(reshape(X{aa}, DimX{aa}(1), prod(DimX{aa}(2:end))).^2));
    SSY(aa)=sum(sum(reshape(Y{aa}, DimY{aa}(1), prod(DimY{aa}(2:end))).^2));
    ssx{aa}=[];
    ssy{aa}=[];
end

for a=1:nLV
   iter = 0;
   %randn('state',sum(clock*100));
   [PCcoeff, PCvec] = pcaKM(Y{1}, 1);
   Tu(:,a) = Y{1}*PCvec;
   Tt(:,a) = randn(size(Tt(:,a)));
   t_old = Tt(:,a)*100;

%    ff =  figure;
   while (sum((t_old - Tt(:,a)).^2) > tol) && (iter < maxiter)
%      figure(ff); plot(iter, sum((t_old - Tt(:,a)).^2), 'o'); hold on;
     iter = iter + 1;
     t_old = Tt(:,a);
     T = zeros(DimX{1}(1), nbX);
     for aa=1:nbX
        tempX = reshape(Xres{aa}, DimX{aa}(1), prod(DimX{aa}(2:end)));
        Wb{aa}(:, a) = (tempX'*Tu(:, a) / (Tu(:, a)'*Tu(:, a)));
        Wb{aa}(:, a) = Wb{aa}(:, a)/norm(Wb{aa}(:, a));
        % compute Tb
        Tb{aa}(:, a) = tempX * Wb{aa}(:, a)/(Wb{aa}(:, a)'*Wb{aa}(:, a));
        T(:, aa) = Tb{aa}(:, a);
     end
     Wt(:,a) = T'*Tu(:,a)/(Tu(:,a)'*Tu(:,a));
     Wt(:,a) = Wt(:,a)/norm(Wt(:,a));
     Tt(:,a) = T*Wt(:,a)/(Wt(:,a)'*Wt(:,a));
     U = zeros(DimY{1}(1), nbX);
     for aa=1:nbX
        Qb{aa}(:,a) = Yres{aa}'*Tt(:,a)/(Tt(:,a)'*Tt(:,a));
        % NEW
% % %         Qb{aa}(:,a) = Yres{aa}'*Tb{aa}(:,a)/(Tb{aa}(:,a)'*Tb{aa}(:,a));
        % END NEW
        Ub{aa}(:,a) = Yres{aa}*Qb{aa}(:,a)/(Qb{aa}(:,a)'*Qb{aa}(:,a));
        U(:, aa) = Ub{aa}(:, a);
     end
     Wu(:,a) = U'*Tt(:,a)/(Tt(:,a)'*Tt(:,a));
     Wu(:,a) = Wu(:,a)/norm(Wu(:,a));
     Tu(:,a) = U*Wu(:,a)/(Wu(:,a)'*Wu(:,a));
   end
   fprintf('number of iterations: %g \n',iter);

   if iter >= maxiter
      s = ['WARNING: maximum number of iterations (' num2str(maxiter) ') reached before convergence'];
      disp(s)
   end


    for aa = 1:nbX        
        % compure B
        B{aa}(1:a, a) = pinv(Tt(:,1:a).'*Tt(:, 1:a))*Tt(:, 1:a).'*Ub{aa}(:, a);
%         B{aa}(1:a, a) = pinv(Tb{aa}(:, 1:a).'*Tb{aa}(:, 1:a))*Tb{aa}(:, 1:a).'*Ub{aa}(:, a);

        % update X and Y
%         Xres{aa} = Xres{aa} - reshape(Tb{aa}(:, a) * kron(W2b{aa}(:, a), W1b{aa}(:, a)).', DimX);
%         Yres{aa} = Y{aa} - Tb{aa}(:, 1:a)*B{aa}(1:a, 1:a)*Qb{aa}(:, 1:a).';
      
        Pb{aa}(:,a) = Xres{aa}'*Tt(:,a)/(Tt(:,a)'*Tt(:,a));
        % NEW
        Wb_reproj{aa}(:,a)  = Xres{aa}'*Tt(:, a) / (Tt(:, a)'*Tt(:, a));
        Wb_reproj{aa}(:,a) = Wb_reproj{aa}(:, a)/norm(Wb_reproj{aa}(:, a));
       % END NEW
        Xres{aa} = Xres{aa} - Tt(:, a) * Pb{aa}(:, a).';
        
% % %         w = lars(normalize(Xres{aa}),center(Tt(:,a)), 'lasso');
% % %         Wridge{aa}(:, a) = w(max(find(sum(isnan(w), 2) == 0)), :).';
% % %         Wridge{aa}(:,a) = Wridge{aa}(:, a)/norm(Wridge{aa}(:, a));

%          xmodel = zeros(size(X{aa}));
%          for z = 1:a
%              xmodel = xmodel+reshape(Tt(:, z) * kron(W2b{aa}(:, z), W1b{aa}(:, z)).', DimX);
%          end
%          Xres{aa} = X{aa}-xmodel;
         

         Yres{aa} = Y{aa} - Tt(:, 1:a)*Qb{aa}(:,1:a).';
         ssx{aa} = [ssx{aa}; ...
             1 - sum(sum(reshape(Xres{aa}, DimX{aa}(1), prod(DimX{aa}(2:end))).^2)) / SSX(aa) ]; 
         ssy{aa} = [ssy{aa}; ...
             1 - sum(sum(reshape(Yres{aa}, DimY{aa}(1), prod(DimY{aa}(2:end))).^2)) / SSY(aa) ];
    end

end

