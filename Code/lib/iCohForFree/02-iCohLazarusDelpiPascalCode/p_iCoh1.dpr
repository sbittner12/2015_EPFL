program p_iCoh1;

{$IFDEF FPC}
  {$MODE Delphi}
{$ENDIF}

uses
{$IFnDEF FPC}
{$ELSE}
  Interfaces,
{$ENDIF}
  Forms,
  u_iCoh1 in 'u_iCoh1.pas' {Form1},
  utils1 in 'utils1.pas',
  uMatToolsDynDouble33 in 'uMatToolsDynDouble33.pas',
  u_iCoh3Math in 'u_iCoh3Math.pas',
  uDFTdouble in 'uDFTdouble.PAS';

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TForm1, Form1);
  Application.Run;
end.
