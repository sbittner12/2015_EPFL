object Form1: TForm1
  Left = 0
  Top = 0
  Caption = 'isolated effective coherence'
  ClientHeight = 307
  ClientWidth = 676
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  Position = poScreenCenter
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object FileTxt: TLabel
    Left = 91
    Top = 20
    Width = 19
    Height = 13
    Caption = 'file?'
    OnClick = FileButtClick
  end
  object Label1: TLabel
    Left = 310
    Top = 97
    Width = 89
    Height = 13
    Caption = '(=0 for automatic)'
  end
  object Label2: TLabel
    Left = 460
    Top = 97
    Width = 151
    Height = 13
    Caption = '(for output results, power of 2)'
  end
  object Label3: TLabel
    Left = 460
    Top = 116
    Width = 106
    Height = 13
    Caption = '(from zero to Nyquist)'
  end
  object FileButt: TButton
    Left = 10
    Top = 15
    Width = 75
    Height = 25
    Caption = 'Select file'
    TabOrder = 0
    OnClick = FileButtClick
  end
  object Ptxt: TLabeledEdit
    Left = 10
    Top = 70
    Width = 121
    Height = 21
    EditLabel.Width = 104
    EditLabel.Height = 13
    EditLabel.Caption = 'Number of time series'
    TabOrder = 1
    Text = '5'
  end
  object Qtxt: TLabeledEdit
    Left = 160
    Top = 70
    Width = 121
    Height = 21
    EditLabel.Width = 43
    EditLabel.Height = 13
    EditLabel.Caption = 'AR order'
    TabOrder = 2
    Text = '3'
  end
  object NTtxt: TLabeledEdit
    Left = 310
    Top = 70
    Width = 121
    Height = 21
    EditLabel.Width = 114
    EditLabel.Height = 13
    EditLabel.Caption = 'Number of time samples'
    TabOrder = 3
    Text = '0'
  end
  object NFreqTxt: TLabeledEdit
    Left = 460
    Top = 70
    Width = 121
    Height = 21
    EditLabel.Width = 119
    EditLabel.Height = 13
    EditLabel.Caption = 'Number of discrete freqs'
    TabOrder = 4
    Text = '128'
  end
  object GoButt: TButton
    Left = 10
    Top = 127
    Width = 75
    Height = 25
    Caption = 'Go'
    TabOrder = 5
    OnClick = GoButtClick
  end
  object OpenDialog1: TOpenDialog
    Filter = 'text|*.txt|all|*.*'
    Left = 529
    Top = 189
  end
end
