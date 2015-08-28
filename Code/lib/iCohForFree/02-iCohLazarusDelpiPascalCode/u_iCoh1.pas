
unit u_iCoh1;

{$IFDEF FPC}
  {$MODE Delphi}
{$ENDIF}

interface

uses
{$IFnDEF FPC}
  Windows,
{$ELSE}
  LCLIntf, LCLType, LMessages,
{$ENDIF}
  Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, ExtCtrls, FileUtil, math, utils1, uMatToolsDynDouble33,
  u_iCoh3Math;

type
  TForm1 = class(TForm)
    FileButt: TButton;
    FileTxt: TLabel;
    OpenDialog1: TOpenDialog;
    Ptxt: TLabeledEdit;
    Qtxt: TLabeledEdit;
    NTtxt: TLabeledEdit;
    Label1: TLabel;
    NFreqTxt: TLabeledEdit;
    Label2: TLabel;
    GoButt: TButton;
    Label3: TLabel;
    procedure FileButtClick(Sender: TObject);
    procedure GoButtClick(Sender: TObject);
    procedure FormCreate(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
    p,q,nt,nfreq:integer;
    fn:string;
    x:matrix;
    function init1:boolean;
  end;

var
  Form1: TForm1;

implementation

{$IFnDEF FPC}
  {$R *.dfm}
{$ELSE}
  {$R *.lfm}
{$ENDIF}

procedure TForm1.FormCreate(Sender: TObject);
begin
  fn:='';
end;

procedure TForm1.FileButtClick(Sender: TObject);
begin
  if opendialog1.Execute then begin
    filetxt.Caption:=opendialog1.FileName;
    fn:=opendialog1.FileName;
  end;
end;

function TForm1.init1:boolean;
var i,j:integer; s:string;
begin
  init1:=false;

  val(ptxt.Text,p,i);
  if i<>0 then begin
    showmessage('error in number of time series');
    exit;
  end;
  if p<2 then begin
    showmessage('error in number of time series < 2');
    exit;
  end;

  val(qtxt.Text,q,i);
  if i<>0 then begin
    showmessage('error in AR order');
    exit;
  end;
  if q<1 then begin
    showmessage('error in AR order <1');
    exit;
  end;

  val(nttxt.Text,nt,i);
  if i<>0 then begin
    showmessage('error in number of time samples');
    exit;
  end;
  if nt<0 then begin
    showmessage('error in number of time samples <0');
    exit;
  end;

  val(nfreqtxt.Text,nfreq,i);
  if i<>0 then begin
    showmessage('error in number of discrete frequencies');
    exit;
  end;
  if nfreq<8 then begin
    showmessage('error in number of discrete frequencies <8');
    exit;
  end;
  i:=round(ln(nfreq)/ln(2));
  j:=round(power(2,i));
  if j<>nfreq then begin
    s:='error in number of discrete frequencies, must be power of two:'+chr(13);
    s:=s+'8 or 16 or 32 or 64 or 128 or 256 or ...';
    showmessage(s);
    exit;
  end;

  if fn='' then begin
    showmessage('select an input file');
    exit;
  end;
  if not FileExistsUTF8(fn) { *Converted from FileExists* } then begin
    showmessage('error in file: can''t find');
    exit;
  end;

  if nt=0 then begin
    nt:=TxtNNs(fn) div p;
    nttxt.Text:=inttostr(nt);
    application.ProcessMessages;
  end;

  nfreq:=2*nfreq;

  readtxt1(nt,p,fn,x);

  init1:=true;
end;

procedure TForm1.GoButtClick(Sender: TObject);
begin
  if not init1 then exit;
  iCoh_Do1(fn,x,q,nfreq);
  showmessage('finished');
end;

end.
