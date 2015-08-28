function [bb, dim, acmat] = iCOH(Data, p, nfreqs)

%% Data: size nt = number of timee point * nv =  number of channel
%% p: order of the function
nt = size(Data,1); %% number of samples
nv = size(Data,2); %% number of variables

%% procedure AutoCov3
dim = nv*(p+1);
acmat = zeros(nv*(p+1),nv*(p+1)); %% create1(dim,dim,acmat); %% fillchar1(acmat,0);
v = zeros(1,nv*(p+1)); %%  create1(dim,v);

for t = 1:nv
    Data_sub(:,t) = Data(:,t) - mean(Data(:,t)); %% subtractcolumnmeans(x) -> Columns of x will have zero mean
end
  
for t = p+1:nt 
    for tt= t:-1:t-p
        k = t-tt;
        for i = 1:nv 
            v(1,i+k*nv) = Data_sub(tt,i);
        end
 
    end
    for i = 1:dim
        for j = i:dim
            acmat(i,j)= acmat(i,j)+v(i)*v(j); %% updates(v,acmat);
        end
    end
    
end


for i =1:dim 
    acmat(i,i) = acmat(i,i)/nt;
    for j = i+1:dim
      acmat(i,j) = acmat(i,j)/nt;
      acmat(j,i) = acmat(i,j);
    end
end

%% acmat is the matrix of the autoregressive model


%% procedure AutoCovMatToPartial
for i = nv+1:dim
    
    %% ok:=SWEEP(i,-1,acmat);
    
    if acmat(i,i) > 0
        T = 1.0/acmat(i,i);
        acmat(i,i)=0.0;
        for jj =1:i-1
            A(jj)=acmat(jj,i);
            acmat(jj,i) =0.0;
        end
        for jj = i+1:size(acmat,1)
            A(jj)=acmat(i,jj);
            acmat(i,jj)=0.0;
        end
        A(i)=-1;
        for ii = 1:size(acmat,1)
            for jj = ii:size(acmat,1)
                acmat(ii,jj)= acmat(ii,jj)-A(ii)*A(jj)*T;
            end
        end
    end
   
end
for i =1:dim
    for j = i+1:dim
      acmat(j,i) = acmat(i,j);
    end
end

%% procedure AutoCovMatToFunc3

for k =1:p+1
    for i = 1:nv
      for j = 1:nv
        acfunc{k}(i,j) = acmat(i,j+(k-1)*nv);
      end
    end
end

%% procedure Autoco2FFT
LN2I=1.4426950408889634073599246810019;
for i = 1:nv
    for j = 1:nv 
        vr = zeros(nfreqs,1); %fillchar1(vr,0);
        if i == j 
            vr(1) =1;
        end
    
        for t = 2 : p+1
          vr(t) = - acfunc{t}(i,j);
        end
        vi = zeros(nfreqs,1); %fillchar1(vi,0);
        vi2 = fft(vr); %% questo è da controllare perchè non sono sicura!
        
%         %% fft
%         n = length(vr);
%         nv2 =n/2;
%         JJ = 1;
%         for II = 1:n-1
%             if II<JJ 
%               TR = vr(JJ);
%               vr(JJ)=vr(II);
%               vr(II)=TR;
%               TI = vi(JJ);
%               vi(JJ) = vi(II);
%               vi(II) = TI;
%             end
%             k=nv2;
%             
%             while k<JJ
%                 JJ=JJ-k;
%                 k = k/2;
%             end
%             JJ= JJ+k;
%         end
%         m = round(LN2I*log(n));
%         LE = 1;
%         for L = 1:m
%             LE=LE*2;
%             LE1=LE/2;
%             UR=1;
%             UI=0;
%             ANG=pi/LE1;
%             WR=cos(ANG);
%             WI=-sin(ANG);
%             for JJ=1:LE1
%                 II=JJ;
%                 while II<=n 
%                     IP=II+LE1;
%                     TR=vr(IP)*UR-vi(IP)*UI;
%                     TI=vr(IP)*UI+vi(IP)*UR;
%                     vr(IP)=vr(II)-TR;
%                     vi(IP)=vi(II)-TI;
%                     vr(II)=vr(II)+TR;
%                     vi(II)=vi(II)+TI;
%                     II=II+LE;
%                 end
%                   TR=UR*WR-UI*WI;
%                   UI=UR*WI+UI*WR;
%                   UR=TR;
%             end
%         end
    
  
        
        for t = 2:nfreqs 
            ccrssj2iRe{t}(i,j) = real(vi2(t));
            ccrssj2iIm{t}(i,j) = imag(vi2(t));
        end
    end
end

ccrssj2iRe{1} = zeros(nv,nv); %fillchar1(ccrssj2iRe[1],0);
ccrssj2iIm{1} = zeros(nv,nv); %fillchar1(ccrssj2iIm[1],0);

%% procedure iCoh_0
ccrss2{1} = zeros(nv,nv);
for t = 2:nfreqs
    for j = 1:nv
        ccrss2{t}(j,j) =0;
    end
    for j = 1:nv
        rj = [];
        rj = (ccrssj2iIm{t}(j,j).^2 + ccrssj2iRe{t}(j,j).^2)/acfunc{1}(j,j);
       
        for i = 1: nv 
            if i<j || i > j
                
                ri = (ccrssj2iIm{t}(i,j).^2 + ccrssj2iRe{t}(i,j).^2)/acfunc{1}(i,i);
                ccrss2{t}(i,j) = ri/(ri+rj); %% ccrss2 dovrebbe essere iCoh_0
            end
        end
    end
    
end

%% procedure WriteTxt1x


for k = 1:size(ccrss2,2)/2
    
            
    bb(k,:) = reshape(ccrss2{k},1,nv*nv);

end

        