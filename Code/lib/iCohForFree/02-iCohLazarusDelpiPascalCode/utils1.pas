
unit utils1;

{$IFDEF FPC}
  {$MODE Delphi}
{$ENDIF}

interface

uses sysutils;

function TxtNNs(fn:string):integer;

implementation




function TxtNNs(fn:string):integer;
label 1,2;
var f:textfile; tt,t:integer; r:real;
begin
  {$I-}
  assignfile(f,fn);
  filemode:=0;
  reset(f);
  tt:=0;
  for t:=1 to maxint do begin
    if eof(f) then begin
      2:;
      tt:=t;
      goto 1;
    end;
    read(f,r);
    if ioresult<>0 then
      goto 2
    else begin
      if r>0 then;
    end;
  end;
  1:;
  TxtNNs:=tt-1;
  closefile(f);
  {$I+}
end;




end.

