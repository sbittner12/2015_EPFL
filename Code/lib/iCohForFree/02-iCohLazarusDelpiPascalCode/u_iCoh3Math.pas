
{$R-}

unit u_iCoh3Math;

{$IFDEF FPC}
  {$MODE Delphi}
{$ENDIF}

interface

uses uMatToolsDynDouble33, uDFTdouble;

procedure iCoh_Do1(savefn:string;
                   var x:matrix; // nt*nv; DC set to zero inside
                   p,nfreqs:integer // p:AR order; nfreqs pow of 2; nfreqs>=p+1
                   );

implementation

uses dialogs, sysutils, math;

var acmat:matrix; acfunc,ccrssj2iRe,ccrssj2iIm,ccrss2:tmatmat;




procedure AutoCov3(var x:matrix; // nt*nv; DC set to zero inside
                   var acmat:matrix; // [nv*(p+1)]*[nv*(p+1)]; created here if cc.nr=0
                   nt,nv,p:integer // p:AR order
                   );
// auto-covariance matrix order: t, t-1, t-2,...,t-p
// acmat[1..nv] is present; acmat[nv+1..2nv] is one step in past; and so on...
var i,tt,t,dim,k,j:integer; v:vector;
begin
  dim:=nv*(p+1);
  if acmat.nr=0 then
    create1(dim,dim,acmat);
  create1(dim,v);
  fillchar1(acmat,0);
  subtractcolumnmeans(x);
  for t:=p+1 to nt do begin
    for tt:=t downto t-p do begin
      k:=t-tt;
      for i:=1 to nv do
        v.v[i+k*nv]:=x.m[tt,i];
    end;
    updates(v,acmat);
  end;
  destroy1(v);
  for i:=1 to dim do begin
    acmat.m[i,i]:=acmat.m[i,i]/nt;
    for j:=i+1 to dim do begin
      acmat.m[i,j]:=acmat.m[i,j]/nt;
      acmat.m[j,i]:=acmat.m[i,j];
    end;
  end;
end;




procedure AutoCovMatToFunc3(nv,p:integer; var acmat:matrix; var acfunc:tmatmat);
// creates acfunc inside, size (p+1)*sqr(nv); p is AR order
var k,i,j:integer;
begin
  for k:=1 to (p+1) do
    create1(nv,nv,acfunc[k]);

  for k:=1 to p+1 do
    for i:=1 to nv do
      for j:= 1 to nv do
        acfunc[k].m[i,j]:=acmat.m[i,j+(k-1)*nv];

end;




procedure AutoCovMatToPartial(var acmat:matrix; // [nv*(p+1)]*[nv*(p+1)] input/output
                              nv,p:integer // p:AR order
                              );
// auto-covariance matrix order: t, t-1, t-2,...,t-p
// acmat[1..nv] is present; acmat[nv+1..2nv] is one step in past; and so on...
// sweeps from nv+1...nv*(p+1), thus leaves AR coefficients un upper right rectangle
var i,j:integer; ok:boolean;
begin
  for i:=(nv+1) to (nv*(p+1)) do begin
    ok:=SWEEP(i,-1,acmat);
    if not ok then begin
      showmessage('Failed partial inverse sweep# '+inttostr(i));
      halt;
    end;
  end;
  for i:=1 to nv*(p+1) do
    for j:=i+1 to nv*(p+1) do
      acmat.m[j,i]:=acmat.m[i,j];
end;




procedure WriteTxt1x(fn:string; nummat:integer; var x:tmatmat);
var txt:textfile; k,i,j:integer; s:string;
begin
  {$I+}
  try
    assignfile(txt,fn);
    filemode:=1;
    rewrite(txt);

    for i:=1 to x[1].nr do
      for j:=1 to x[1].nc do begin
        s:='['+inttostr(i)+'>'+inttostr(j)+']';
        write(txt,s:16);
      end;
    writeln(txt);

    for k:=1 to nummat div 2 do begin
      for j:=1 to x[k].nr do
        for i:=1 to x[k].nc do
          write(txt,' ',x[k].m[i,j]:15);
      writeln(txt);
    end;
  finally
    closeFile(txt);
  end;
end;




procedure Autoco2FFT(nv,nfreqs,p:integer; var acfunc,ccrssj2iRe,ccrssj2iIm:tmatmat);
// input auto-covar funtion[1..p+1], output causal cross-spectra ccrss created
// inside [1..nfreqs] power of 2
// BaccalaSameshima version
var i,j,t:integer; vr,vi:vector;
begin
  for i:=1 to nfreqs do begin // power of 2, nfreqs>=(p+1)
    create1(nv,nv,ccrssj2iRe[i]);
    create1(nv,nv,ccrssj2iIm[i]);
  end;
  create1(nfreqs,vr);
  create1(nfreqs,vi);

  for i:=1 to nv do
    for j:=1 to nv do begin
      fillchar1(vr,0);
      if i=j then
        vr.v[1]:=1;
      for t:=2 to p+1 do
        vr.v[t]:=-acfunc[t].m[i,j];
      fillchar1(vi,0);
      fft(vr,vi);
      for t:=2 to nfreqs do begin
        ccrssj2iRe[t].m[i,j]:=vr.v[t];
        ccrssj2iIm[t].m[i,j]:=vi.v[t];
      end;
    end;

  fillchar1(ccrssj2iRe[1],0);
  fillchar1(ccrssj2iIm[1],0);

  destroy1(vr);
  destroy1(vi);

end;




procedure iCoh_0(nv,nfreqs:integer; var cov:matrix; var ccrssj2iRe,ccrssj2iIm,ccrss2:tmatmat);
/// non scaled, just one global factor;
var i,j,t:integer; ri,rj:extended;
begin
  for i:=1 to nfreqs do
    create1(nv,nv,ccrss2[i]);
  fillchar1(ccrss2[1],0);
  for t:=2 to nfreqs do begin
    for j:=1 to nv do
      ccrss2[t].m[j,j]:=0;
    for j:=1 to nv do begin
      rj:=(sqr(ccrssj2iRe[t].m[j,j])+sqr(ccrssj2iIm[t].m[j,j]))/cov.m[j,j];
      for i:=1 to nv do if i<>j then begin
        ri:=(sqr(ccrssj2iRe[t].m[i,j])+sqr(ccrssj2iIm[t].m[i,j]))/cov.m[i,i];;
        ccrss2[t].m[i,j]:=ri/(ri+rj);
      end;
    end;
  end;
end;




procedure gPDC(nv,nfreqs:integer; var cov:matrix; var ccrssj2iRe,ccrssj2iIm,ccrss2:tmatmat);
var i,j,t:integer; r,r1:extended;
begin
  fillchar1(ccrss2[1],0);
  for t:=2 to nfreqs do begin
    for j:=1 to nv do begin
      r:=0;
      for i:=1 to nv do begin
        r1:=(sqr(ccrssj2iRe[t].m[i,j])+sqr(ccrssj2iIm[t].m[i,j]))/cov.m[i,i];
        ccrss2[t].m[i,j]:=r1;
        r:=r+r1;
      end;
      if r<=0 then r:=1;
      for i:=1 to nv do begin
        ccrss2[t].m[i,j]:=ccrss2[t].m[i,j]/r;
      end;
    end
  end;
  for t:=2 to nfreqs do
    for j:=1 to nv do
      ccrss2[t].m[j,j]:=0;
end;




procedure iCoh_Do1(savefn:string;
                   var x:matrix; // nt*nv; DC set to zero inside
                   p,nfreqs:integer // p:AR order; nfreqs pow of 2; nfreqs>=p+1
                   );
var nt,nv:integer; fn:string;
begin
  if nfreqs<(p+1) then exit;
  acmat.nr:=0;
  nt:=x.nr;
  nv:=x.nc;
  AutoCov3(x,acmat,nt,nv,p); // later destroy acmat
  AutoCovMatToPartial(acmat,nv,p);
  AutoCovMatToFunc3(nv,p,acmat,acfunc); // later destroy p+1, acfunc
  fn:=changefileext(savefn,'-Noisevar-ARcoeff.txt');
  WriteTxt1(fn,p+1,acfunc);
  Autoco2FFT(nv,nfreqs,p,acfunc,ccrssj2iRe,ccrssj2iIm); // later destroy nfreqs ccrss
  iCoh_0(nv,nfreqs,acfunc[1],ccrssj2iRe,ccrssj2iIm,ccrss2);
  fn:=changefileext(savefn,'-iCoh.txt');
  WriteTxt1x(fn,nfreqs,ccrss2);

  gPDC(nv,nfreqs,acfunc[1],ccrssj2iRe,ccrssj2iIm,ccrss2);
  fn:=changefileext(savefn,'-gPDC.txt');
  WriteTxt1x(fn,nfreqs,ccrss2);

  destroy3(nfreqs,ccrss2);
  destroy3(p+1,acfunc);
  destroy3(nfreqs,ccrssj2iRe);
  destroy3(nfreqs,ccrssj2iIm);
  destroy1(acmat);
end;




end.

