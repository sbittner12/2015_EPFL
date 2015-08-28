function [Tb,W1b,W2b,Wt,Tt,Ub,Qb,Wu,Tu, ssx, ssy, W1b_temp, W2b_temp, Tu_orig] = MBtriPLS2_microstates(X,Y,microstates,nLV,tol)
% function [Tb,Pb,Wb,Wt,Tt,Ub,Qb,Wu,Tu,B,ssq,Rbo,Rbv,Ry,Lbo,Lbv,Lto] =
% MBtriPLS2_cmmnW1b(X,Y,nLV,tol)
% in  :
% X  (1 x number of X-blocks) cell array where each element is an X-data-block 
%   of dimension t (time) x w (frequency) x d (channel) 
% Y  (1 x number of Y-blocks) cell array where each element is an Y-data-block 
%   of dimension t (time) x imf
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
% B (X-variables x Y-variables.max(nLV)) regression vectors
% ssq (max(nLV) x 1+number of X-blocks+1+number of Y-blocks) cumulative sum of squares, [Xt X-blocks Yt Y-blocks]
% Rbo (objects x max(nLV)) X-block object residuals
% Rbv (all variables x max(nLV)) X-block variable residuals
% Ry (objects x Y-variables.max(nLV)) Y-block residuals
% Lbo (objects x max(nLV)) X-block object leverages
% Lbv (all variables x max(nLV)) X-block variable leverages
% Lto (objects x max(nLV)) X-block super score object leverages

if nargin == 0
   help MBspls
   return;
elseif nargin == 3
   tol = 1e-12;
end

maxiter = 1000;
nbX = length(X);
DimX = size(X{1});
DimY = size(Y{1});

Xres = X;
Yres = Y;
Tb = cell(1, nbX); 
% Pb = cell(1, nbX); 
W1b = cell(1, nbX); 
W2b = microstates;
W1b_temp = cell(1, nbX); 
W2b_temp = cell(1, nbX); 
for b=1:nbX
    Tb{b} = zeros(DimX(1), nLV);
    W1b{b} = zeros(DimX(2), nLV);
end
Wt = zeros(nbX,nLV);
Tt = zeros(DimX(1), nLV);
Ub = cell(1, nbX); 
Qb = cell(1, nbX); 
for b = 1:nbX
    Ub{b} = zeros(DimY(1), nLV);
    Qb{b} = zeros(DimY(2), nLV);
end
Wu = zeros(nbX,nLV);
Tu = zeros(DimX(1),nLV);
Tu_orig = zeros(DimX(1), nLV);
B = cell(1, nbX);
for b = 1:nbX
    B{b} = zeros(nLV,nLV);
end


SSX = zeros(1, nbX);
SSY = zeros(1, nbX);
ssx = cell(1, nbX);
ssy = cell(1, nbX);
for aa = 1:nbX
    SSX(aa)=sum(sum(reshape(X{aa}, DimX(1), prod(DimX(2:end))).^2));
    SSY(aa)=sum(sum(reshape(Y{aa}, DimY(1), prod(DimY(2:end))).^2));
    ssx{aa}=[];
    ssy{aa}=[];
end

for a=1:nLV
   iter = 0;
   %randn('state',sum(clock*100));
   Tu(:,a) = Y{1}(:,1);
   Tu_orig(:,a) = Tu(:,a);
   Tt(:,a) = randn(size(Tt(:,a)));
   t_old = Tt(:,a)*100;

%    ff =  figure;
   while (sum((t_old - Tt(:,a)).^2) > tol) && (iter < maxiter)
%      figure(ff); plot(iter, sum((t_old - Tt(:,a)).^2), 'o'); hold on;
     iter = iter + 1;
     t_old = Tt(:,a);
     T = zeros(DimX(1), nbX);
     Wb_reshaped = cell(1, nbX);
     for aa=1:nbX
        tempX = reshape(Xres{aa}, DimX(1), prod(DimX(2:end)));
        Wb = tempX'*Tu(:,a);
        Wb_reshaped{aa} = reshape(Wb,DimX(2),DimX(3));
     end
     Wb_all = cell2mat(Wb_reshaped);
     [u, s, v] = svd(Wb_all);
     for aa = 1:nbX
        W1b{aa}(:, a) = (u(:, 1));
%         W2b{aa}(:, a) = v( (aa-1)*DimX(3)+1:aa*DimX(3), 1);
        W1b{aa}(:, a) = W1b{aa}(:, a)/norm(W1b{aa}(:, a));
%         W2b{aa}(:, a) = W2b{aa}(:, a)/norm(W2b{aa}(:, a));
        % compute Tb
        Tb{aa}(:, a) = tempX * kron(W2b{aa}(:, a), W1b{aa}(:, a));
        T(:, aa) = Tb{aa}(:, a);
     end
     Wt(:,a) = T'*Tu(:,a)/(Tu(:,a)'*Tu(:,a));
     Wt(:,a) = Wt(:,a)/norm(Wt(:,a));
     Tt(:,a) = T*Wt(:,a)/(Wt(:,a)'*Wt(:,a));
     U = zeros(DimY(1), nbX);
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

   if iter == maxiter
      s = ['WARNING: maximum number of iterations (' num2str(maxiter) ') reached before convergence'];
      disp(s)
   end

    % NEW: May 13, 2010
    % Recompute W2b (spatial patterns)
% %     Wb_reshaped = cell(1, nbX);
% %     for aa=1:nbX
% %         tempX = reshape(Xres{aa}, DimX(1), prod(DimX(2:end)));
% %         Wb = tempX'*Tt(:,a);
% %         Wb_reshaped{aa} = reshape(Wb,DimX(2),DimX(3));
% %      end
% %      Wb_all = cell2mat(Wb_reshaped);
% %      [u, s, v] = svd(Wb_all);
% %      for aa = 1:nbX
% %         W1b_temp{aa}(:, a) = (u(:, 1));
% %         W2b_temp{aa}(:, a) = v( (aa-1)*DimX(3)+1:aa*DimX(3), 1);
% %         W1b_temp{aa}(:, a) = W1b_temp{aa}(:, a)/norm(W1b_temp{aa}(:, a));
% %         W2b_temp{aa}(:, a) = W2b_temp{aa}(:, a)/norm(W2b_temp{aa}(:, a));
% %      end
    for aa = 1:nbX
        comp = Tb{aa}(:, a) * W1b{aa}(:, a).';
        W2b_temp{aa}(:, a) = ...
            reshape(Xres{aa}, prod(DimX(1:2)), DimX(3)).' ...
            * reshape(comp, [], 1);
        W2b_temp{aa}(:, a) = W2b_temp{aa}(:, a)/norm(W2b_temp{aa}(:, a));
    end

    for aa = 1:nbX         
        
        % compure B
        B{aa}(1:a, a) = pinv(Tt(:,1:a).'*Tt(:, 1:a))*Tt(:, 1:a).'*Ub{aa}(:, a);
%         B{aa}(1:a, a) = pinv(Tb{aa}(:, 1:a).'*Tb{aa}(:, 1:a))*Tb{aa}(:, 1:a).'*Ub{aa}(:, a);

        % update X and Y
%         Xres{aa} = Xres{aa} - reshape(Tb{aa}(:, a) * kron(W2b{aa}(:, a), W1b{aa}(:, a)).', DimX);
%         Yres{aa} = Y{aa} - Tb{aa}(:, 1:a)*B{aa}(1:a, 1:a)*Qb{aa}(:, 1:a).';

         Xres{aa} = Xres{aa} - reshape(Tt(:, a) * kron(W2b{aa}(:, a), W1b{aa}(:, a)).', DimX);

%          xmodel = zeros(size(X{aa}));
%          for z = 1:a
%              xmodel = xmodel+reshape(Tt(:, z) * kron(W2b{aa}(:, z), W1b{aa}(:, z)).', DimX);
%          end
%          Xres{aa} = X{aa}-xmodel;
         

         Yres{aa} = Y{aa} - Tt(:, 1:a)*B{aa}(1:a, 1:a)*Qb{aa}(:,1:a).';
         ssx{aa} = [ssx{aa}; ...
             1 - sum(sum(reshape(Xres{aa}, DimX(1), prod(DimX(2:end))).^2)) / SSX(aa) ]; 
         ssy{aa} = [ssy{aa}; ...
             1 - sum(sum(reshape(Yres{aa}, DimY(1), prod(DimY(2:end))).^2)) / SSY(aa) ];
    end
    


%    if a > 1
%       index = (a-1)*p+1:a*p;
%       B(:,index) = B(:,index-p);
%    end
%    index = (a-1)*p+1:a*p;
%    B(:,index) = B(:,index) + Wb(:,a)*inv(Pb(:,a)'*Wb(:,a))*Qb(:,a)';
      
%    ssq(a,1) = (ssqXY(1) - sum(sum(X.^2)))/ssqXY(1);
%    ssq(a,nbX+2) = (ssqXY(nbX+2) - sum(sum(Y.^2)))/ssqXY(nbX+2);
%    for aa=1:nbX
%       rowi = Xin(aa,1):Xin(aa,2);
%       coli = (a-1)*nbX+aa;
%       ssq(a,aa+1) = (ssqXY(aa+1) - sum(sum(X(:,rowi).^2)))/ssqXY(aa+1);
%       Rbo(:,coli) = sqrt(sum(X(:,rowi).^2,2));
% %       index = aa:nbX:(a-1)*nbX+aa;
% %       Lbo(:,coli) = diag(Tb(:,index)*pinv(Tb(:,index)'*Tb(:,index))*Tb(:,index)');
%    end
%    for aa=1:nbY
%       rowi = Yin(aa,1):Yin(aa,2);
%       coli = (a-1)*nbY+aa;
%       ssq(a,aa+2+nbX) = (ssqXY(aa+2+nbX) - sum(sum(Y(:,rowi).^2)))/ssqXY(aa+2+nbX);
%       Ry(:,coli) = sqrt(sum(Y(:,rowi).^2,2));
%    end
%    Rbv(:,a) = sqrt(sum(X.^2,1))';
%    Lbv(:,a) = diag(Pb(:,1:a)*Pb(:,1:a)');
%    Lto(:,a) = diag(Tt(:,1:a)*pinv(Tt(:,1:a)'*Tt(:,1:a))*Tt(:,1:a)');
end

