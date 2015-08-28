
unit uMatToolsDynDouble33;

{$IFDEF FPC}
  {$MODE Delphi}
{$ENDIF}

interface

uses sysutils, FileUtil;

{$R-}

(* GENERAL NOTE:
You must turn range checking off "{$R-}" if you want to access the elements
of vectors and matrices explicitly in your program!
*)


type real=double; // change here if you like {single double extended}

var CheckingOn:boolean=false;
    Ok_NoError:boolean;
    kk1:integer=1;
    kk2:integer=2;
    kk3:integer=3;
    kk4:integer=4;
{Set CheckingOn:=true to force checking of correct vector/matrix dimensions.
Then you will have to manually check the value of Ok_NoError by yourself, since
the program will not stop by error! Error checking is left up to you if you set
CheckingOn:=true!!}


type TVec3=array[1..3] of real;
     TInternalVector=array[1..1] of real;
     pvector=^TInternalVector;
     vector=record
       n:integer;
       v:pvector;
     end;
     TInternalMatrix=array[1..1] of pvector;
     pmatrix=^TInternalMatrix;
     matrix=record
       nr,nc:integer;
       m:pmatrix;
     end;
     TInternalIntVector=array[1..1] of integer;
     pintvector=^TInternalIntVector;
     intvector=record
       n:integer;
       v:pintvector;
     end;
     TInternalIntMatrix=array[1..1] of pintvector;
     pintmatrix=^TInternalIntMatrix;
     intmatrix=record
       nr,nc:integer;
       m:pintmatrix;
     end;

     // for an array of matrices
     TInternalMatrix3=array[1..1] of matrix;
     pmatrix3=^TInternalMatrix3;
     matrix3=record
       nm,nr,nc:integer;
       m3:pmatrix3;
     end;

     // Explicit array of matrices
     TMatMat=array[1..8196] of matrix;

     // Explicit matrix of matrices
     TMatOfMat=array[1..1024,1..512] of matrix;

// For explicit array of matrices
procedure Create3(nm,nr,nc:integer; var m:TMatMat); // also sets all to 0
procedure Destroy3(nm:integer; var m:TMatMat);

// For explicit matrix of matrices
procedure Create5(OuterNR,OuterNC,nr,nc:integer; var m:TMatOfMat); // also sets all to 0
procedure Destroy5(OuterNR,OuterNC:integer; var m:TMatOfMat);

// for an array of matrices
procedure Create2(nm,nr,nc:integer; var m:matrix3);
procedure Destroy2(var m:matrix3);

{ It's your absolute responsibility to create and destroy variables, using
these procedures}
procedure Create1(n:integer; var v:vector); overload;
procedure Create1(nr,nc:integer; var m:matrix); overload;
procedure Create1(n:integer; var v:intvector); overload;
procedure Create1(nr,nc:integer; var m:intmatrix); overload;
procedure Destroy1(var v:vector); overload;
procedure Destroy1(var m:matrix); overload;
procedure Destroy1(var v:intvector); overload;
procedure Destroy1(var m:intmatrix); overload;

{ GENERAL NOTES on accessing elements of vectors/matrices.
1. If v is a pvector, then v[i] is of type real. E.g., if a and b reals,
then the following expressions are valid: a:=b+v[5] and v[10]:=v[8]-a.
There's no need for writing v^[5] under delphi 7.
2. If a and b are pmatrices, then a[i,j] and b[i,j] are of type real;
no need for writing a^[i]^[j]]^ etc. a[1,5]-b[6,3] is valid.
3. if a is a pmatrix, then a^[i] is of type pvector. Here we do need the "^" symbol!
}


{NOTE: Never equate two vectors or two matrices. This must be done element
by element or using the "Equate" procedures.}
procedure Equate1(a,b:matrix); overload; // a=b
procedure Equate1(a,b:vector); overload; // a=b
procedure Equate1(a,b:intvector); overload; // a=b
procedure Equate1(a,b:intmatrix); overload; // a=b
procedure Equate3(nm:integer; a,b:tmatmat); overload; // a=b
{-----------------------------------------------}


{NOTE: Never use fillchar on vectors or matrices. This must be done element
by element or using the "Fill" procedures.}
procedure FillChar1(var a:vector; val:real); overload; // a=val
procedure FillChar1(var a:matrix; val:real); overload; // a=val
procedure FillChar1(var a:intvector; val:integer); overload; // a=val
procedure FillChar1(var a:intmatrix; val:integer); overload; // a=val
{-----------------------------------------------}


function VxV(n:integer; var a,b:pvector):extended; overload; // plain scalar product

function VxV(var a,b:vector):extended; overload; // plain scalar product

procedure MxV(var a:vector; var b:matrix; var c:vector); overload; // a=bc; creates a if a.n=0

procedure MxV(var a:vector; var b:matrix; var c:pvector);  overload; // a=bc; creates a if a.n=0

procedure VxM(var r,v:vector; var m:matrix); overload; // r=v*m; creates r if r.n=0

procedure VxM(var r:vector; var v:pvector; var m:matrix);  overload; // r=v*m; creates r if r.n=0

procedure MxM(var a,b,c:matrix); // a=bc; creates a if a.nr=0

procedure MplusM(var a,b,c:matrix); // a=b+c; creates a if a.nr=0

procedure MxMt(var a,m:matrix); overload; // a=m*m'; creates a if a.nr=0

procedure MxMt(var a,b,c:matrix); overload; // a=bc'; creates a if a.nr=0

procedure MtxM(var a,m:matrix); overload; // a=m'*m; creates a if a.nr=0

procedure MtxM(var a,b,c:matrix); overload; // a=b'c; creates a if a.nr=0

function Trace(var m:matrix):real;

function UnitTrace(var m:matrix):real;
// returns tr and makes m=m/tr (forces matrix to unit trace)

function SWEEP(INDEX:integer;
               SWPDIR:integer; { -1 FOR FORWARD SWEEP, +1 OTHERWISE }
               VAR M:MATRIX;
               {temp work space below}
               var a:vector // if a.n=0 then creates and destroys internally
               ):boolean; overload;
{Only upper triangle used. If you sweep for index=1..order, with swpdir=-1,
you get "-inverse" in the upper triangle.}

function SWEEP(INDEX:integer;
               SWPDIR:integer; { -1 FOR FORWARD SWEEP, +1 OTHERWISE }
               VAR M:MATRIX
               ):boolean; overload;
{Only upper triangle used. If you sweep for index=1..order, with swpdir=-1,
you get "-inverse" in the upper triangle. Creates & destroys an internal vector.}

procedure RowNorms(var ns:vector; var x:matrix);
{ns will contain norms of rows of x. Creates ns if ns.n=0}

procedure ColumnNorms(var ns:vector; var x:matrix);
{ns will contain norms of columns of x. Creates ns if ns.n=0}

procedure UpdateS(var v:vector; var s:matrix); overload;
{Initialize s=0. Performs update s=s+vv' only for upper triangle.}

procedure UpdateS(f:extended; var v:vector;var s:matrix); overload;
{Initialize s=0. Performs update s=s+f*vv' only for upper triangle.}

procedure UpdateS(var v:pvector; var s:matrix); overload;
{Initialize s=0. Performs update s=s+vv' only for upper triangle.}

procedure UpdateS(f:extended; var v:pvector; var s:matrix); overload;
{Initialize s=0. Performs update s=s+f*vv' only for upper triangle.}

procedure UpdateS(var x,y:pvector; var s:matrix); overload;
{Initialize s=0. Performs update s=s+x*y' all matrix elements.}

procedure UpdateS(f:extended; var x,y:pvector; var s:matrix); overload;
{Initialize s=0. Performs update s=s+x*y' all matrix elements.}

procedure RowMeans(var ns:vector; var x:matrix);
{ns will contain means of rows of x. Creates ns if ns.n=0}

procedure SubtractMean(var x:vector); overload;
{x will have zero mean}

procedure SubtractRowMeans(var ns:vector; var x:matrix); overload;
{Rows of x will have zero mean. ns will contain means of rows of x.
Creates ns if ns.n=0}

procedure SubtractRowMeans(var x:matrix); overload;
{Rows of x will have zero mean.}

procedure SubtractRowMeans(nr,nc:integer; var x:pmatrix); overload;
{Rows of x will have zero mean.}

procedure SubtractColumnMeans(var x:matrix); overload;
{Columns of x will have zero mean.}

procedure SubtractColumnMeans(nr,nc:integer; var x:pmatrix); overload;
{Columns of x will have zero mean.}

procedure StandardizeColumns(var x:matrix); overload;
{Columns of x will have zero mean and unit variance.}

procedure StandardizeColumns(var x:matrix; var sd:vector); overload;
{Columns of x will have zero mean and unit variance. outputs sd, creates if sd.n=0}

procedure NormalizeV(var x:vector);
{x will have unit norm}

function NormalizeVf(var x:vector):extended; overload;
{x will have unit norm, returns norm}

procedure NormalizeV2(var x:vector);
{x will have average sum of squares = 1}

procedure NormalizeColumns(var x:matrix);
{Columns of x will have unit norms}

procedure NormalizeColumns1(var x:matrix);
{For non-negative x, columns of x will have average sum = 1}

procedure NormalizeColumns2(var x:matrix);
{Columns of x will have squared norms = nr (number of rows)}

procedure StandardizeRows(var x:matrix);
{Rows of x will have zero mean and unit variance.}

procedure NormalizeRows(var x:matrix);
{Rows of x will have unit norms}

procedure NormalizeRows1(var x:matrix);
{For non-negative x, rows of x will have average sum = 1}

procedure NormalizeMatrix(var x:matrix);
{Matrix will have average sum of squares = 1}

procedure NormalizeMatrix1(var x:matrix);
{For non-negative x, matrix will have average sum = 1}

procedure SubtractGrandMean(var x:matrix);
{x will have zero grand mean.}

function ReadTxt1(nr,nc:integer; fn:string; var x:matrix):boolean; overload;
{Creates and reads in text matrix. Last column can be any string, but will be
ignored. Footer allowed and ignored.}

function ReadTxt1(n:integer; fn:string; var v:vector):boolean; overload;
{Creates and reads in text vector}

function ReadTxt1(fn:string; var x:matrix):boolean; overload;
{Reads in text matrix (previously created). Last column can be any string,
but will be ignored. Footer allowed and ignored.}

function ReadTxtT1(nr,nc:integer; fn:string; var xt:matrix):boolean;
{Creates and reads in text matrix transposed, i.e., in file nr*nc, but placed
into xt nc*nr. Last column can be any string, but will be ignored. Footer
allowed and ignored.}

function ReadReal1(nr,nc:integer; fn:string; var x:matrix):boolean; overload;
{Creates and reads in real matrix}

function ReadReal1(fn:string; var x:matrix):boolean; overload;
{reads in real matrix}

function ReadByte1(nr,nc:integer; fn:string; var x:matrix):boolean;
{Creates and reads in byte matrix, e.g., raw grey image.}

function WriteTxt1(fn:string; var x:matrix):boolean; overload;
{Writes text matrix.}

function WriteTxt1append(fn:string; var x:matrix):boolean; overload;
{Writes text matrix appending to existing file. If file don't exist, then creates file.}

function WriteTxt1(fn:string; nummat:integer; var x:tmatmat):boolean; overload;
{Writes text many matrices.}

function WriteTxt1(fn:string; var x:vector):boolean; overload;
{Writes text vector.}

function WriteSingle1(fn:string; var x:vector):boolean; overload;
{Writes single vector.}

function WriteSingle2(fn:string; var x:vector):boolean;
{Writes single vector, appending at end of file, if already exists.}

function WriteSingle3(fn:string; var x:matrix):boolean; overload;
{Writes single matrix.}

function WriteTxtT1(fn:string; var x:matrix):boolean;
{Writes transposed text matrix.}

function WriteReal1(fn:string; var x:matrix):boolean; overload;
{Writes real matrix.}

function WriteReal1(fn:string; var x:vector):boolean; overload;
{Writes real vector}

function WriteRealT1(fn:string; var x:matrix):boolean;
{Writes transposed real matrix.}

function WriteByte1(fn:string; var x:matrix):boolean; overload;
{Writes byte matrix, e.g., raw grey image.}

function WriteByteT1(fn:string; var x:matrix):boolean;
{Writes transposed byte matrix, e.g., raw grey image.}

procedure AveRef(var v:vector);

function ReadFromSingle1(nr,nc:integer; fn:string; var x:matrix):boolean; overload;
{Creates and reads in real matrix from file of single}

function ReadFromSingle1(n:integer; fn:string; var v:vector):boolean; overload;
{Creates and reads in real vector from file of single}

implementation

uses math;

var eps1:extended=0.0;
    eps2:extended=0.0;
    precis:extended absolute eps1;
    eta:extended absolute eps2;




procedure Create1(n:integer; var v:vector); overload;
begin
  if n=0 then exit;
  v.n:=n;
  getmem(v.v,n*sizeof(real));
end;

procedure Destroy1(var v:vector); overload;
begin
  if v.n=0 then exit;
  freemem(v.v,v.n*sizeof(real));
  v.n:=0;
end;

procedure Create1(n:integer; var v:intvector); overload;
begin
  if n=0 then exit;
  v.n:=n;
  getmem(v.v,n*sizeof(integer));
end;

procedure Destroy1(var v:intvector); overload;
begin
  if v.n=0 then exit;
  freemem(v.v,v.n*sizeof(integer));
  v.n:=0;
end;

procedure Create1(nr,nc:integer; var m:matrix); overload;
var i:integer;
begin
  if nr=0 then exit;
  m.nr:=nr;
  m.nc:=nc;
  getmem(m.m,nr*sizeof(pointer));
  for i:=1 to nr do
    getmem(m.m[i],nc*sizeof(real))
end;

procedure Destroy1(var m:matrix); overload;
var i:integer;
begin
  if m.nr=0 then exit;
  for i:=m.nr downto 1 do
    freemem(m.m[i],m.nc*sizeof(real));
  freemem(m.m,m.nr*sizeof(pointer));
  m.nr:=0;
  m.nc:=0;
end;

procedure Create1(nr,nc:integer; var m:intmatrix); overload;
var i:integer;
begin
  if nr=0 then exit;
  m.nr:=nr;
  m.nc:=nc;
  getmem(m.m,nr*sizeof(pointer));
  for i:=1 to nr do
    getmem(m.m[i],nc*sizeof(integer))
end;

procedure Destroy1(var m:intmatrix); overload;
var i:integer;
begin
  if m.nr=0 then exit;
  for i:=m.nr downto 1 do
    freemem(m.m[i],m.nc*sizeof(integer));
  freemem(m.m,m.nr*sizeof(pointer));
  m.nr:=0;
  m.nc:=0;
end;




procedure Check1(var a,b:matrix); overload;
begin
  Ok_NoError:=(a.nr=b.nr) and (a.nc=b.nc);
end;

procedure Check1(var a,b:intmatrix); overload;
begin
  Ok_NoError:=(a.nr=b.nr) and (a.nc=b.nc);
end;

procedure Check1(var a,b:vector); overload;
begin
  Ok_NoError:=(a.n=b.n)
end;

procedure Check1(var a,b:intvector); overload;
begin
  Ok_NoError:=(a.n=b.n)
end;

procedure Check1(var a:vector; var b:matrix; var c:vector); overload;
begin
  Ok_NoError:=(b.nc=c.n);
  if not Ok_NoError then exit;
  if a.n>0 then
    Ok_NoError:=(a.n=b.nr);
end;

procedure Check1(var r,v:vector; var m:matrix); overload;
begin
  Ok_NoError:=(m.nr=v.n);
  if not Ok_NoError then exit;
  if r.n>0 then
    Ok_NoError:=(r.n=m.nc);
end;

procedure Check1(var a,b,c:matrix); overload;
begin
  Ok_NoError:=(b.nc=c.nr);
  if not Ok_NoError then exit;
  if a.nr>0 then
    Ok_NoError:=(a.nr=b.nr) and (a.nc=c.nc);
end;




procedure Equate1(a,b:vector); overload; // a=b
var i:integer;
begin
  if CheckingOn then begin
    Check1(a,b);
    if not Ok_NoError then exit;
  end;
  for i:=1 to a.n do
    a.v[i]:=b.v[i];
end;

procedure Equate1(a,b:matrix); overload; // a=b
var i,j:integer;
begin
  if CheckingOn then begin
    Check1(a,b);
    if not Ok_NoError then exit;
  end;
  for i:=1 to a.nr do
    for j:=1 to a.nc do
      a.m[i,j]:=b.m[i,j];
end;

procedure Equate3(nm:integer; a,b:tmatmat); overload; // a=b
var k:integer;
begin
  for k:=1 to nm do
    equate1(a[k],b[k]);
end;

procedure Equate1(a,b:intvector); overload; // a=b
var i:integer;
begin
  if CheckingOn then begin
    Check1(a,b);
    if not Ok_NoError then exit;
  end;
  for i:=1 to a.n do
    a.v[i]:=b.v[i];
end;

procedure Equate1(a,b:intmatrix); overload; // a=b
var i,j:integer;
begin
  if CheckingOn then begin
    Check1(a,b);
    if not Ok_NoError then exit;
  end;
  for i:=1 to a.nr do
    for j:=1 to a.nc do
      a.m[i,j]:=b.m[i,j];
end;




procedure FillChar1(var a:matrix; val:real); overload; // a=val
var i,j:integer;
begin
  for i:=1 to a.nr do
    for j:=1 to a.nc do
      a.m[i,j]:=val;
end;

procedure FillChar1(var a:intmatrix; val:integer); overload; // a=val
var i,j:integer;
begin
  for i:=1 to a.nr do
    for j:=1 to a.nc do
      a.m[i,j]:=val;
end;

procedure FillChar1(var a:vector; val:real); overload; // a=val
var i:integer;
begin
  for i:=1 to a.n do
    a.v[i]:=val;
end;

procedure FillChar1(var a:intvector; val:integer); overload; // a=val
var i:integer;
begin
  for i:=1 to a.n do
    a.v[i]:=val;
end;




function VxV(n:integer; var a,b:pvector):extended; overload;
var t:extended; i:integer;
begin
  t:=0.0;
  for i:=1 to n do
    t:=t+a[i]*b[i];
  VxV:=t
end;




function VxV(var a,b:vector):extended; overload;
var t:extended; i:integer;
begin
  VxV:=0;
  if CheckingOn then begin
    Check1(a,b);
    if not Ok_NoError then exit;
  end;
  t:=0.0;
  for i:=1 to a.n do
    t:=t+a.v[i]*b.v[i];
  VxV:=t;
end;




procedure MxV(var a:vector; var b:matrix; var c:vector); overload; // a=bc; creates a if a.n=0
var i:integer;
begin
  if CheckingOn then begin
    Check1(a,b,c);
    if not Ok_NoError then exit;
  end;
  if a.n=0 then
    Create1(b.nr,a);
  for i:=1 to a.n do
    a.v[i]:=VxV(b.nc,b.m^[i],c.v);
end;




procedure MxV(var a:vector; var b:matrix; var c:pvector);  overload; // a=bc; creates a if a.n=0
var i:integer;
begin
  if a.n=0 then
    Create1(b.nr,a);
  for i:=1 to a.n do
    a.v[i]:=VxV(b.nc,b.m^[i],c);
end;




procedure VxM(var r,v:vector; var m:matrix); overload; // r=v*m; creates r if r.n=0
VAR I,J:INTEGER; T:extended;
BEGIN { VECxMAT }
  if CheckingOn then begin
    Check1(r,v,m);
    if not Ok_NoError then exit;
  end;
  if r.n=0 then
    Create1(m.nc,r);
  FOR J:=1 TO m.nc DO BEGIN
    T:=0.0;
    FOR I:=1 TO m.nr DO
      T:=T+V.v[I]*M.m[I,J];
    R.v[J]:=T
  END
END; { VECxMAT }




procedure VxM(var r:vector; var v:pvector; var m:matrix); overload; // r=v*m; creates r if r.n=0
VAR I,J:INTEGER; T:extended;
BEGIN { VECxMAT }
  if r.n=0 then
    Create1(m.nc,r);
  FOR J:=1 TO m.nc DO BEGIN
    T:=0.0;
    FOR I:=1 TO m.nr DO
      T:=T+V[I]*M.m[I,J];
    R.v[J]:=T
  END
END; { VECxMAT }




procedure MxM(var a,b,c:matrix); // a=bc; creates a if a.nr=0
var i,j,k:integer; t:extended;
begin
  if CheckingOn then begin
    Check1(a,b,c);
    if not Ok_NoError then exit;
  end;
  if a.nr=0 then
    Create1(b.nr,c.nc,a);
  for i:=1 to a.nr do
    for j:=1 to a.nc do begin
      t:=0;
      for k:=1 to b.nc do
        t:=t+b.m[i,k]*c.m[k,j];
      a.m[i,j]:=t
    end;
end;




procedure MplusM(var a,b,c:matrix); // a=b+c; creates a if a.nr=0
var i,j:integer; t:extended;
begin
  if a.nr=0 then
    Create1(b.nr,b.nc,a);
  for i:=1 to a.nr do
    for j:=1 to a.nc do begin
      t:=b.m[i,j]+c.m[i,j];
      a.m[i,j]:=t;
    end;
end;




procedure MxMt(var a,m:matrix); overload; // a=m*m'; creates a if a.nr=0
var i,j:integer;
begin

  if CheckingOn then
    if a.nr>0 then begin
      Ok_NoError:=(a.nr=m.nr) and ((a.nc=m.nr));
      if not Ok_NoError then exit;
    end;
  if a.nr=0 then
    Create1(m.nr,m.nr,a);

  for i:=1 to a.nr do
    for j:=i to a.nc do
      a.m[i,j]:=vxv(m.nc,m.m^[i],m.m^[j]);
  for i:=1 to a.nr do
    for j:=i+1 to a.nc do
      a.m[j,i]:=a.m[i,j];
end;




procedure MxMt(var a,b,c:matrix); overload; // a=bc'; creates a if a.nr=0
var i,j:integer;
begin
  if CheckingOn then begin
    Ok_NoError:=(b.nc=c.nc);
    if not Ok_NoError then exit;
    if a.nr>0 then begin
      Ok_NoError:=(a.nr=b.nr) and ((a.nc=c.nr));
      if not Ok_NoError then exit;
    end;
  end;
  if a.nr=0 then
    Create1(b.nr,c.nr,a);

  for i:=1 to a.nr do
    for j:=1 to a.nc do
      a.m[i,j]:=vxv(b.nc,b.m^[i],c.m^[j]);
end;




procedure MtxM(var a,m:matrix); overload; // a=m'*m; creates a if a.nr=0
var i,j,k:integer; t:extended;
begin
  if CheckingOn then
    if a.nr>0 then begin
      Ok_NoError:=(a.nr=m.nc) and ((a.nc=m.nc));
      if not Ok_NoError then exit;
    end;
  if a.nr=0 then
    Create1(m.nc,m.nc,a);

  for i:=1 to a.nr do
    for j:=i to a.nc do begin
      t:=0;
      for k:=1 to m.nr do
        t:=t+m.m[k,i]*m.m[k,j];
      a.m[i,j]:=t
    end;
  for i:=1 to a.nr do
    for j:=i+1 to a.nc do
      a.m[j,i]:=a.m[i,j];
end;




procedure MtxM(var a,b,c:matrix); overload; // a=b'c; creates a if a.nr=0
var i,j,k:integer; t:extended;
begin
  if CheckingOn then begin
    Ok_NoError:=(b.nr=c.nr);
    if not Ok_NoError then exit;
    if a.nr>0 then begin
      Ok_NoError:=(a.nr=b.nc) and ((a.nc=c.nc));
      if not Ok_NoError then exit;
    end;
  end;
  if a.nr=0 then
    Create1(b.nc,c.nc,a);

  for i:=1 to a.nr do
    for j:=1 to a.nc do begin
      t:=0;
      for k:=1 to b.nr do
        t:=t+b.m[k,i]*c.m[k,j];
      a.m[i,j]:=t
    end;
end;




function SWEEP1(ORDER,INDEX:integer;
                SWPDIR:integer; { -1 FOR FORWARD SWEEP, +1 OTHERWISE }
                VAR M:pMATRIX;
                {temp work space below}
                var a:pvector // order*1
                ):boolean;
{ Only upper triangle used. If you sweep for index=1..order, with swpdir=-1,
you get "-inverse" in the upper triangle.}
VAR I,J:integer; T:extended;
BEGIN { SWEEP1 }
  sweep1:=false;

  if abs(m[index,index])<eps2 then
    exit;

  T:=1.0/M[INDEX,INDEX];
  M[INDEX,INDEX]:=0.0;
  FOR I:=1 TO INDEX-1 DO BEGIN
    A[I]:=M[I,INDEX];
    M[I,INDEX]:=0.0
  END;
  FOR I:=INDEX+1 TO ORDER DO BEGIN
    A[I]:=M[INDEX,I];
    M[INDEX,I]:=0.0
  END;
  A[INDEX]:=SWPDIR;
  FOR I:=1 TO ORDER DO
    FOR J:=I TO ORDER DO
      M[I,J]:=M[I,J]-A[I]*A[J]*T;
  sweep1:=true;

END; { SWEEP1 }




function SWEEP(INDEX:integer;
               SWPDIR:integer; { -1 FOR FORWARD SWEEP, +1 OTHERWISE }
               VAR M:MATRIX;
               {temp work space below}
               var a:vector // if a.n=0 then creates and destroys internally
               ):boolean; overload;
{Only upper triangle used. If you sweep for index=1..order, with swpdir=-1,
you get "-inverse" in the upper triangle.}
VAR cr:boolean;
BEGIN { SWEEP }
  cr:=(a.n=0);
  if cr then
    create1(m.nr,a);
  if CheckingOn then
    if a.n<>m.nr then begin
      Ok_NoError:=false;
      sweep:=false;
      if cr then
        destroy1(a);
      exit;
    end;
  sweep:=SWEEP1(m.nr,INDEX,SWPDIR,m.m,a.v);
  if cr then
    destroy1(a);
END; { SWEEP }




function SWEEP(INDEX:integer;
               SWPDIR:integer; { -1 FOR FORWARD SWEEP, +1 OTHERWISE }
               VAR M:MATRIX
               ):boolean; overload;
{Only upper triangle used. If you sweep for index=1..order, with swpdir=-1,
you get "-inverse" in the upper triangle. Creates & destroys an internal vector.}
VAR a:vector;
BEGIN { SWEEP }
  create1(m.nr,a);
  sweep:=SWEEP1(m.nr,INDEX,SWPDIR,m.m,a.v);
  destroy1(a);
END; { SWEEP }




procedure TRANSP(VAR Mt,M:MATRIX); // mt is m transposed; creates mt if mt.nr=0
VAR I,J:integer;
BEGIN { TRANSP }
  if mt.nr=0 then
    create1(m.nc,m.nr,mt);
  if CheckingOn then
    if (mt.nr<>m.nc) or (mt.nc<>m.nr) then begin
      Ok_NoError:=false;
      exit;
    end;
  FOR I:=1 TO m.nr DO
    FOR J:=1 TO m.nc DO
      mt.m[J,I]:=M.m[I,J]
END; { TRANSP }




procedure RowNorms(var ns:vector; var x:matrix);
{ns will contain norms of rows of x. Creates ns if ns.n=0}
var i:integer;
begin
  if CheckingOn then
    if ns.n>0 then begin
      Ok_NoError:=(x.nr=ns.n);
      if not Ok_NoError then exit;
    end;
  if ns.n=0 then
    Create1(x.nr,ns);
  for i:=1 to ns.n do
    ns.v[i]:=sqrt(vxv(x.nc,x.m^[i],x.m^[i]));
end;




procedure ColumnNorms(var ns:vector; var x:matrix);
{ns will contain norms of columns of x. Creates ns if ns.n=0}
var i,j:integer; r1:extended;
begin
  if CheckingOn then
    if ns.n>0 then begin
      Ok_NoError:=(x.nc=ns.n);
      if not Ok_NoError then exit;
    end;
  if ns.n=0 then
    Create1(x.nc,ns);
  for i:=1 to ns.n do begin
    r1:=0;
    for j:=1 to x.nr do
      r1:=r1+sqr(x.m[j,i]);
    ns.v[i]:=sqrt(r1);
  end;
end;




procedure UpdateS(var v:vector;
                  var s:matrix
                  ); overload;
{Initialize s=0. Performs update s=s+vv' only for upper triangle.}
var i,j:integer;
begin { UpdateS }
  if CheckingOn then begin
    Ok_NoError:=(v.n=s.nr) and (v.n=s.nc);
    if not Ok_NoError then exit;
  end;
  for i:=1 to v.n do
    for j:=i to v.n do
      s.m[i,j]:=s.m[i,j]+v.v[i]*v.v[j];
end; { UpdateS }




procedure UpdateS(f:extended;
                  var v:vector;
                  var s:matrix
                  ); overload;
{Initialize s=0. Performs update s=s+f*vv' only for upper triangle.}
var i,j:integer;
begin { UpdateS }
  if CheckingOn then begin
    Ok_NoError:=(v.n=s.nr) and (v.n=s.nc);
    if not Ok_NoError then exit;
  end;
  for i:=1 to v.n do
    for j:=i to v.n do
      s.m[i,j]:=s.m[i,j]+f*v.v[i]*v.v[j];
end; { UpdateS }




procedure UpdateS(var v:pvector;
                  var s:matrix
                  ); overload;
{Initialize s=0. Performs update s=s+vv' only for upper triangle.}
var i,j:integer;
begin { UpdateS }
  if CheckingOn then begin
    Ok_NoError:=(s.nr=s.nc);
    if not Ok_NoError then exit;
  end;
  for i:=1 to s.nr do
    for j:=i to s.nc do
      s.m[i,j]:=s.m[i,j]+v[i]*v[j];
end; { UpdateS }




procedure UpdateS(f:extended;
                  var v:pvector;
                  var s:matrix
                  ); overload;
{Initialize s=0. Performs update s=s+f*vv' only for upper triangle.}
var i,j:integer;
begin { UpdateS }
  if CheckingOn then begin
    Ok_NoError:=(s.nr=s.nc);
    if not Ok_NoError then exit;
  end;
  for i:=1 to s.nr do
    for j:=i to s.nc do
      s.m[i,j]:=s.m[i,j]+f*v[i]*v[j];
end; { UpdateS }




procedure UpdateS(var x,y:pvector;
                  var s:matrix
                  ); overload;
{Initialize s=0. Performs update s=s+x*y' all matrix elements.}
var i,j:integer;
begin { UpdateS }
  for i:=1 to s.nr do
    for j:=1 to s.nc do
      s.m[i,j]:=s.m[i,j]+x[i]*y[j];
end; { UpdateS }




procedure UpdateS(f:extended;
                  var x,y:pvector;
                  var s:matrix
                  ); overload;
{Initialize s=0. Performs update s=s+x*y' all matrix elements.}
var i,j:integer;
begin { UpdateS }
  for i:=1 to s.nr do
    for j:=1 to s.nc do
      s.m[i,j]:=s.m[i,j]+f*x[i]*y[j];
end; { UpdateS }




procedure RowMeans(var ns:vector; var x:matrix);
{ns will contain means of rows of x. Creates ns if ns.n=0}
var i,j:integer; r1:extended;
begin
  if CheckingOn then
    if ns.n>0 then begin
      Ok_NoError:=(x.nr=ns.n);
      if not Ok_NoError then exit;
    end;
  if ns.n=0 then
    Create1(x.nr,ns);
  for i:=1 to x.nr do begin
    r1:=0;
    for j:=1 to x.nc do
      r1:=r1+x.m[i,j];
    ns.v[i]:=r1/x.nc;
  end;
end;




procedure SubtractRowMeans(var ns:vector; var x:matrix); overload;
{Rows of x will have zero mean. ns will contain means of rows of x.
Creates ns if ns.n=0}
var i,j:integer; r1:extended;
begin
  if CheckingOn then
    if ns.n>0 then begin
      Ok_NoError:=(x.nr=ns.n);
      if not Ok_NoError then exit;
    end;
  if ns.n=0 then
    Create1(x.nr,ns);
  for i:=1 to x.nr do begin
    r1:=0;
    for j:=1 to x.nc do
      r1:=r1+x.m[i,j];
    r1:=r1/x.nc;
    ns.v[i]:=r1;
    for j:=1 to x.nc do
      x.m[i,j]:=x.m[i,j]-r1;
  end;
end;




procedure SubtractRowMeans(var x:matrix); overload;
{Rows of x will have zero mean.}
var i,j:integer; r1:extended;
begin
  for i:=1 to x.nr do begin
    r1:=0;
    for j:=1 to x.nc do
      r1:=r1+x.m[i,j];
    r1:=r1/x.nc;
    for j:=1 to x.nc do
      x.m[i,j]:=x.m[i,j]-r1;
  end;
end;




procedure SubtractRowMeans(nr,nc:integer; var x:pmatrix); overload;
{Rows of x will have zero mean.}
var i,j:integer; r1:extended;
begin
  for i:=1 to nr do begin
    r1:=0;
    for j:=1 to nc do
      r1:=r1+x[i,j];
    r1:=r1/nc;
    for j:=1 to nc do
      x[i,j]:=x[i,j]-r1;
  end;
end;




procedure SubtractColumnMeans(var x:matrix); overload;
{Columns of x will have zero mean.}
var i,j:integer; r1:extended;
begin
  for j:=1 to x.nc do begin
    r1:=0;
    for i:=1 to x.nr do
      r1:=r1+x.m[i,j];
    r1:=r1/x.nr;
    for i:=1 to x.nr do
      x.m[i,j]:=x.m[i,j]-r1;
  end;
end;




procedure SubtractColumnMeans(nr,nc:integer; var x:pmatrix); overload;
{Columns of x will have zero mean.}
var i,j:integer; r1:extended;
begin
  for j:=1 to nc do begin
    r1:=0;
    for i:=1 to nr do
      r1:=r1+x[i,j];
    r1:=r1/nr;
    for i:=1 to nr do
      x[i,j]:=x[i,j]-r1;
  end;
end;




procedure StandardizeColumns(var x:matrix); overload;
{Columns of x will have zero mean and unit variance.}
var i,j:integer; r1,r2:extended;
begin
  for j:=1 to x.nc do begin
    r1:=0;
    r2:=0;
    for i:=1 to x.nr do begin
      r1:=r1+x.m[i,j];
      r2:=r2+sqr(x.m[i,j]);
    end;
    r1:=r1/x.nr;
    r2:=sqrt(r2/x.nr-sqr(r1));
    for i:=1 to x.nr do
      x.m[i,j]:=(x.m[i,j]-r1)/r2;
  end;
end;




procedure StandardizeColumns(var x:matrix; var sd:vector); overload;
{Columns of x will have zero mean and unit variance. outputs sd, creates if sd.n=0}
var i,j:integer; r1,r2:extended;
begin
  if sd.n=0 then
    create1(x.nc,sd);
  for j:=1 to x.nc do begin
    r1:=0;
    r2:=0;
    for i:=1 to x.nr do begin
      r1:=r1+x.m[i,j];
      r2:=r2+sqr(x.m[i,j]);
    end;
    r1:=r1/x.nr;
    r2:=sqrt(r2/x.nr-sqr(r1));
    sd.v[j]:=r2;
    for i:=1 to x.nr do
      x.m[i,j]:=(x.m[i,j]-r1)/r2;
  end;
end;




procedure NormalizeColumns(var x:matrix);
{Columns of x will have unit norms}
var i,j:integer; r2:extended;
begin
  for j:=1 to x.nc do begin
    r2:=0;
    for i:=1 to x.nr do
      r2:=r2+sqr(x.m[i,j]);
    if r2<eps2 then
      r2:=1;
    r2:=sqrt(r2);
    for i:=1 to x.nr do
      x.m[i,j]:=x.m[i,j]/r2;
  end;
end;




procedure NormalizeColumns2(var x:matrix);
{Columns of x will have squared norms = nr (number of rows)}
var i,j:integer; r2:extended;
begin
  for j:=1 to x.nc do begin
    r2:=0;
    for i:=1 to x.nr do
      r2:=r2+sqr(x.m[i,j]);
    if r2<eps2 then
      r2:=1;
    r2:=sqrt(r2/x.nr);
    for i:=1 to x.nr do
      x.m[i,j]:=x.m[i,j]/r2;
  end;
end;




procedure NormalizeColumns1(var x:matrix);
{For non-negative x, columns of x will have average sum = 1}
var i,j:integer; r2:extended;
begin
  for j:=1 to x.nc do begin
    r2:=0;
    for i:=1 to x.nr do
      r2:=r2+x.m[i,j];
    r2:=r2/x.nr;
    if r2<eps2 then
      r2:=1;
    for i:=1 to x.nr do
      x.m[i,j]:=x.m[i,j]/r2;
  end;
end;




procedure StandardizeRows(var x:matrix);
{Rows of x will have zero mean and unit variance.}
var i,j:integer; r1,r2:extended;
begin
  for i:=1 to x.nr do begin
    r1:=0;
    r2:=0;
    for j:=1 to x.nc do begin
      r1:=r1+x.m[i,j];
      r2:=r2+sqr(x.m[i,j]);
    end;
    r1:=r1/x.nc;
    r2:=sqrt(r2/x.nc-sqr(r1));
    for j:=1 to x.nc do
      x.m[i,j]:=(x.m[i,j]-r1)/r2;
  end;
end;




procedure NormalizeRows(var x:matrix);
{Rows of x will have unit norms}
var i,j:integer; r2:extended;
begin
  for i:=1 to x.nr do begin
    r2:=0;
    for j:=1 to x.nc do
      r2:=r2+sqr(x.m[i,j]);
    if r2<eps2 then
      r2:=1;
    r2:=sqrt(r2);
    for j:=1 to x.nc do
      x.m[i,j]:=x.m[i,j]/r2;
  end;
end;




procedure NormalizeRows1(var x:matrix);
{For non-negative x, rows of x will have average sum = 1}
var i,j:integer; r2:extended;
begin
  for i:=1 to x.nr do begin
    r2:=0;
    for j:=1 to x.nc do
      r2:=r2+x.m[i,j];
    r2:=r2/x.nc; // average sum of squares
    if r2<eps2 then
      r2:=1;
    for j:=1 to x.nc do
      x.m[i,j]:=x.m[i,j]/r2;
  end;
end;




procedure NormalizeMatrix(var x:matrix);
{Matrix will have average sum of squares = 1}
var i,j:integer; r:extended;
begin
  r:=0;
  for i:=1 to x.nr do
    for j:=1 to x.nc do
      r:=r+sqr(x.m[i,j]);
  r:=r/(x.nr*x.nc); // average sum of squares
  if r<eps2 then
    r:=1;
  r:=sqrt(r);
  for i:=1 to x.nr do
    for j:=1 to x.nc do
      x.m[i,j]:=x.m[i,j]/r;
end;




procedure NormalizeMatrix1(var x:matrix);
{For non-negative x, matrix will have average sum = 1}
var i,j:integer; r:extended;
begin
  r:=0;
  for i:=1 to x.nr do
    for j:=1 to x.nc do
      r:=r+x.m[i,j];
  r:=r/(x.nr*x.nc); // average sum
  if r<eps2 then
    r:=1;
  for i:=1 to x.nr do
    for j:=1 to x.nc do
      x.m[i,j]:=x.m[i,j]/r;
end;




procedure SubtractGrandMean(var x:matrix);
{x will have zero grand mean.}
var i,j:integer; r1:extended;
begin
  r1:=0;
  for i:=1 to x.nr do
    for j:=1 to x.nc do
      r1:=r1+x.m[i,j];
  r1:=r1/(x.nc*x.nr);
  for i:=1 to x.nr do
    for j:=1 to x.nc do
      x.m[i,j]:=x.m[i,j]-r1;
end;




function ReadTxt1(nr,nc:integer; fn:string; var x:matrix):boolean; overload;
{Creates and reads in text matrix. Last column can be any string, but will be
ignored. Footer allowed and ignored.}
var txt:textfile; i,j:integer; s:string;
begin
  {$I+}
  Result:=true;
  create1(nr,nc,x);
  assignfile(txt,fn);
  filemode:=0;
  try
    reset(txt);
    try
      for i:=1 to x.nr do begin
        for j:=1 to x.nc do
          read(txt,x.m[i,j]);
        readln(txt,s);
      end;
    finally
      closeFile(txt);
    end;
  except
    on EInOutError do
      Result:=false;
    on EInvalidOp do
      Result:=false;
  end;
  if not Result then
    destroy1(x);
end;




function ReadTxt1(n:integer; fn:string; var v:vector):boolean; overload;
{Creates and reads in text vector}
var txt:textfile; i:integer;
begin
  {$I+}
  Result:=true;
  create1(n,v);
  assignfile(txt,fn);
  filemode:=0;
  try
    reset(txt);
    try
      for i:=1 to v.n do
        read(txt,v.v[i]);
    finally
      closeFile(txt);
    end;
  except
    on EInOutError do
      Result:=false;
    on EInvalidOp do
      Result:=false;
  end;
  if not Result then
    destroy1(v);
end;




function ReadTxt1(fn:string; var x:matrix):boolean; overload;
{Reads in text matrix (previously created). Last column can be any string,
but will be ignored. Footer allowed and ignored.}
var txt:textfile; i,j:integer; s:string; r:extended;
begin
  {$I+}
  Result:=true;
  assignfile(txt,fn);
  filemode:=0;
  try
    reset(txt);
    try
      for i:=1 to x.nr do begin
        for j:=1 to x.nc do begin
          read(txt,r);
          x.m[i,j]:=r;
        end;
        readln(txt,s);
      end;
    finally
      closeFile(txt);
    end;
  except
    on EInOutError do
      Result:=false;
    on EInvalidOp do
      Result:=false;
  end;
end;




function ReadTxtT1(nr,nc:integer; fn:string; var xt:matrix):boolean;
{Creates and reads in text matrix transposed, i.e., in file nr*nc, but placed
into xt nc*nr. Last column can be any string, but will be ignored. Footer
allowed and ignored.}
var txt:textfile; i,j:integer; s:string;
begin
  {$I+}
  Result:=true;
  create1(nc,nr,xt);
  assignfile(txt,fn);
  filemode:=0;
  try
    reset(txt);
    try
      for i:=1 to nr do begin
        for j:=1 to nc do
          read(txt,xt.m[j,i]);
        readln(txt,s);
      end;
    finally
      closeFile(txt);
    end;
  except
    on EInOutError do
      Result:=false;
    on EInvalidOp do
      Result:=false;
  end;
  if not Result then
    destroy1(xt);
end;




function WriteTxt1(fn:string; var x:matrix):boolean; overload;
{Writes text matrix.}
var txt:textfile; i,j:integer;
begin
  {$I+}
  Result:=true;
  assignfile(txt,fn);
  filemode:=1;
  try
    rewrite(txt);
    try
      for i:=1 to x.nr do begin
        for j:=1 to x.nc do
          write(txt,' ',x.m[i,j]:15);
        writeln(txt);
      end;
    finally
      closeFile(txt);
    end;
  except
    on EInOutError do
      Result:=false;
    on EInvalidOp do
      Result:=false;
  end;
end;




function WriteTxt1append(fn:string; var x:matrix):boolean; overload;
{Writes text matrix appending to existing file. If file don't exist, then creates file.}
var txt:textfile; i,j:integer;
begin
  {$I+}

  if not FileExistsUTF8(fn) { *Converted from FileExists* } then begin

    Result:=true;
    assignfile(txt,fn);
    filemode:=1;
    try
      rewrite(txt);
      try
        for i:=1 to x.nr do begin
          for j:=1 to x.nc do
            write(txt,' ',x.m[i,j]:15);
          writeln(txt);
        end;
      finally
        closeFile(txt);
      end;
    except
      on EInOutError do
        Result:=false;
      on EInvalidOp do
        Result:=false;
    end;

  end else begin

    Result:=true;
    assignfile(txt,fn);
    filemode:=1;
    try
      append(txt);
      try
        for i:=1 to x.nr do begin
          for j:=1 to x.nc do
            write(txt,' ',x.m[i,j]:15);
          writeln(txt);
        end;
      finally
        closeFile(txt);
      end;
    except
      on EInOutError do
        Result:=false;
      on EInvalidOp do
        Result:=false;
    end;

  end;
end;




function WriteTxt1(fn:string; nummat:integer; var x:tmatmat):boolean; overload;
{Writes text many matrices.}
var txt:textfile; k,i,j:integer;
begin
  {$I+}
  Result:=true;
  assignfile(txt,fn);
  filemode:=1;
  try
    rewrite(txt);
    try
      for k:=1 to nummat do begin
        for i:=1 to x[k].nr do begin
          for j:=1 to x[k].nc do
            write(txt,' ',x[k].m[i,j]:15);
          writeln(txt);
        end;
        writeln(txt);
      end;
    finally
      closeFile(txt);
    end;
  except
    on EInOutError do
      Result:=false;
    on EInvalidOp do
      Result:=false;
  end;
end;




function WriteSingle1(fn:string; var x:vector):boolean; overload;
{Writes single vector.}
var txt:file of single; i:integer; s:single;
begin
  {$I+}
  Result:=true;
  assignfile(txt,fn);
  filemode:=1;
  try
    rewrite(txt);
    try
      for i:=1 to x.n do begin
        s:=x.v[i];
        write(txt,s);
      end;
    finally
      closeFile(txt);
    end;
  except
    on EInOutError do
      Result:=false;
    on EInvalidOp do
      Result:=false;
  end;
end;




function WriteSingle3(fn:string; var x:matrix):boolean; overload;
{Writes single matrix.}
var txt:file of single; i,j:integer; s:single;
begin
  {$I+}
  Result:=true;
  assignfile(txt,fn);
  filemode:=1;
  try
    rewrite(txt);
    try
      for i:=1 to x.nr do
        for j:=1 to x.nc do begin
          s:=x.m[i,j];
          write(txt,s);
        end;
    finally
      closeFile(txt);
    end;
  except
    on EInOutError do
      Result:=false;
    on EInvalidOp do
      Result:=false;
  end;
end;




function WriteSingle2(fn:string; var x:vector):boolean;
{Writes single vector, appending at end of file, if already exists.}
var txt:file of single; i:integer; s:single;
begin
  {$I+}
  if not FileExistsUTF8(fn) { *Converted from FileExists* } then begin
    Result:=true;
    assignfile(txt,fn);
    filemode:=1;
    try
      rewrite(txt);
      try
        for i:=1 to x.n do begin
          s:=x.v[i];
          write(txt,s);
        end;
      finally
        closeFile(txt);
      end;
    except
      on EInOutError do
        Result:=false;
      on EInvalidOp do
        Result:=false;
    end;
  end else begin
    Result:=true;
    assignfile(txt,fn);
    filemode:=2;
    try
      reset(txt);
      Seek(txt,FileSize(txt));
      try
        for i:=1 to x.n do begin
          s:=x.v[i];
          write(txt,s);
        end;
      finally
        closeFile(txt);
      end;
    except
      on EInOutError do
        Result:=false;
      on EInvalidOp do
        Result:=false;
    end;
  end;
end;




function WriteTxt1(fn:string; var x:vector):boolean;
{Writes text vector.}
var txt:textfile; i:integer;
begin
  {$I+}
  Result:=true;
  assignfile(txt,fn);
  filemode:=1;
  try
    rewrite(txt);
    try
      for i:=1 to x.n do
        writeln(txt,' ',x.v[i]:15);
    finally
      closeFile(txt);
    end;
  except
    on EInOutError do
      Result:=false;
    on EInvalidOp do
      Result:=false;
  end;
end;




function WriteTxtT1(fn:string; var x:matrix):boolean;
{Writes transposed text matrix.}
var txt:textfile; i,j:integer;
begin
  {$I+}
  Result:=true;
  assignfile(txt,fn);
  filemode:=1;
  try
    rewrite(txt);
    try
      for j:=1 to x.nc do begin
        for i:=1 to x.nr do
          write(txt,' ',x.m[i,j]:15);
        writeln(txt);
      end;
    finally
      closeFile(txt);
    end;
  except
    on EInOutError do
      Result:=false;
    on EInvalidOp do
      Result:=false;
  end;
end;




function ReadReal1(nr,nc:integer; fn:string; var x:matrix):boolean; overload;
{Creates and reads in real matrix}
var txt:file; i:integer;
begin
  {$I+}
  Result:=true;
  create1(nr,nc,x);
  assignfile(txt,fn);
  filemode:=0;
  try
    reset(txt,1);
    try
      for i:=1 to x.nr do
        blockread(txt,x.m^[i]^,sizeof(real)*x.nc);
    finally
      closeFile(txt);
    end;
  except
    on EInOutError do
      Result:=false;
    on EInvalidOp do
      Result:=false;
  end;
  if not Result then
    destroy1(x);
end;




function ReadReal1(fn:string; var x:matrix):boolean; overload;
{reads in real matrix}
var txt:file; i:integer;
begin
  {$I+}
  Result:=true;
  assignfile(txt,fn);
  filemode:=0;
  try
    reset(txt,1);
    try
      for i:=1 to x.nr do
        blockread(txt,x.m^[i]^,sizeof(real)*x.nc);
    finally
      closeFile(txt);
    end;
  except
    on EInOutError do
      Result:=false;
    on EInvalidOp do
      Result:=false;
  end;
end;




function ReadFromSingle1(nr,nc:integer; fn:string; var x:matrix):boolean;
{Creates and reads in real matrix from file of single}
var txt:file; i,j:integer; r:single;
begin
  {$I+}
  Result:=true;
  create1(nr,nc,x);
  assignfile(txt,fn);
  filemode:=0;
  try
    reset(txt,1);
    try
      for i:=1 to x.nr do
        for j:=1 to x.nc do begin
          blockread(txt,r,sizeof(r));
          x.m[i,j]:=r;
        end;
    finally
      closeFile(txt);
    end;
  except
    on EInOutError do
      Result:=false;
    on EInvalidOp do
      Result:=false;
  end;
  if not Result then
    destroy1(x);
end;




function ReadFromSingle1(n:integer; fn:string; var v:vector):boolean; overload;
{Creates and reads in real vector from file of single}
var txt:file; i:integer; r:single;
begin
  {$I+}
  Result:=true;
  create1(n,v);
  assignfile(txt,fn);
  filemode:=0;
  try
    reset(txt,1);
    try
      for i:=1 to v.n do begin
        blockread(txt,r,sizeof(r));
        v.v[i]:=r;
      end;
    finally
      closeFile(txt);
    end;
  except
    on EInOutError do
      Result:=false;
    on EInvalidOp do
      Result:=false;
  end;
  if not Result then
    destroy1(v);
end;




function WriteRealT1(fn:string; var x:matrix):boolean;
{Writes transposed real matrix.}
var txt:file; i,j:integer;
begin
  {$I+}
  Result:=true;
  assignfile(txt,fn);
  filemode:=1;
  try
    rewrite(txt,1);
    try
      for j:=1 to x.nc do
        for i:=1 to x.nr do
          blockwrite(txt,x.m[i,j],sizeof(real));
    finally
      closeFile(txt);
    end;
  except
    on EInOutError do
      Result:=false;
    on EInvalidOp do
      Result:=false;
  end;
end;




function WriteReal1(fn:string; var x:matrix):boolean; overload;
{Writes real matrix}
var txt:file; i:integer;
begin
  {$I+}
  Result:=true;
  assignfile(txt,fn);
  filemode:=1;
  try
    rewrite(txt,1);
    try
      for i:=1 to x.nr do
        blockwrite(txt,x.m^[i]^,x.nc*sizeof(real));
    finally
      closeFile(txt);
    end;
  except
    on EInOutError do
      Result:=false;
    on EInvalidOp do
      Result:=false;
  end;
end;




function WriteReal1(fn:string; var x:vector):boolean; overload;
{Writes real vector}
var txt:file;
begin
  {$I+}
  Result:=true;
  assignfile(txt,fn);
  filemode:=1;
  try
    rewrite(txt,1);
    try
      blockwrite(txt,x.v^,x.n*sizeof(real));
    finally
      closeFile(txt);
    end;
  except
    on EInOutError do
      Result:=false;
    on EInvalidOp do
      Result:=false;
  end;
end;




function ReadByte1(nr,nc:integer; fn:string; var x:matrix):boolean;
{Creates and reads in byte matrix, e.g., raw grey image.}
var txt:file of byte; i,j:integer; b:byte;
begin
  {$I+}
  Result:=true;
  create1(nr,nc,x);
  assignfile(txt,fn);
  filemode:=0;
  try
    reset(txt);
    try
      for i:=1 to x.nr do
        for j:=1 to x.nc do begin
          read(txt,b);
          x.m[i,j]:=b;
        end;
    finally
      closeFile(txt);
    end;
  except
    on EInOutError do
      Result:=false;
    on EInvalidOp do
      Result:=false;
  end;
  if not Result then
    destroy1(x);
end;




function WriteByte1(fn:string; var x:matrix):boolean; overload;
{Writes byte matrix, e.g., raw grey image.}
var txt:file of byte; i,j:integer; b:byte;
begin
  {$I+}
  Result:=true;
  assignfile(txt,fn);
  filemode:=1;
  try
    rewrite(txt);
    try
      for i:=1 to x.nr do
        for j:=1 to x.nc do begin
          if x.m[i,j]>255 then
            b:=255
          else
            if x.m[i,j]<0 then
              b:=0
            else
              b:=round(x.m[i,j]);
          write(txt,b);
        end;
    finally
      closeFile(txt);
    end;
  except
    on EInOutError do
      Result:=false;
    on EInvalidOp do
      Result:=false;
  end;
end;




function WriteByteT1(fn:string; var x:matrix):boolean;
{Writes transposed byte matrix, e.g., raw grey image.}
var txt:file of byte; i,j:integer; b:byte;
begin
  {$I+}
  Result:=true;
  assignfile(txt,fn);
  filemode:=1;
  try
    rewrite(txt);
    try
      for j:=1 to x.nc do
        for i:=1 to x.nr do begin
          if x.m[i,j]>255 then
            b:=255
          else
            if x.m[i,j]<0 then
              b:=0
            else
              b:=round(x.m[i,j]);
          write(txt,b);
        end;
    finally
      closeFile(txt);
    end;
  except
    on EInOutError do
      Result:=false;
    on EInvalidOp do
      Result:=false;
  end;
end;




procedure AveRef(var v:vector);
var r:extended; i:integer;
begin
  r:=0;
  for i:=1 to v.n do
    r:=r+v.v[i];
  r:=r/v.n;
  for i:=1 to v.n do
    v.v[i]:=v.v[i]-r;
end;




procedure AveRefOnColumns(var m:matrix);
{will take each column of m and force average reference on each one}
var r:extended; i,j:integer;
begin
  for j:=1 to m.nc do begin
    r:=0;
    for i:=1 to m.nr do
      r:=r+m.m[i,j];
    r:=r/m.nr;
    for i:=1 to m.nr do
      m.m[i,j]:=m.m[i,j]-r;
  end;
end;




procedure AveRefOnRows(var m:matrix);
{will take each row of m and force average reference on each one}
var r:extended; i,j:integer;
begin
  for i:=1 to m.nr do begin
    r:=0;
    for j:=1 to m.nc do
      r:=r+m.m[i,j];
    r:=r/m.nc;
    for j:=1 to m.nc do
      m.m[i,j]:=m.m[i,j]-r;
  end;
end;




procedure MxD(var a,m:matrix; // a=m*diag(d), a,m are nr*nc
              var d:vector // d is nc
              ); // all pre-created
var i,j:integer;
begin
  for i:=1 to m.nr do
    for j:=1 to m.nc do
      a.m[i,j]:=m.m[i,j]*d.v[j];
end;




procedure Create2(nm,nr,nc:integer; var m:matrix3);
var i:integer;
begin
  m.nm:=nm;
  m.nr:=nr;
  m.nc:=nc;
  getmem(m.m3,nm*sizeof(matrix));
  for i:=1 to nm do
    create1(nr,nc,m.m3[i]);
end;

procedure Destroy2(var m:matrix3);
var i:integer;
begin
  for i:=1 to m.nm do
    destroy1(m.m3[i]);
  freemem(m.m3,m.nm*sizeof(matrix));
  m.nm:=0;
  m.nr:=0;
  m.nc:=0;
end;




procedure Create3(nm,nr,nc:integer; var m:TMatMat);
var k:integer;
begin
  for k:=1 to nm do begin
    Create1(nr,nc,m[k]);
    fillchar1(m[k],0);
  end;
end;




procedure Destroy3(nm:integer; var m:TMatMat);
var k:integer;
begin
  for k:=nm downto 1 do
    destroy1(m[k]);
end;




procedure Create5(OuterNR,OuterNC,nr,nc:integer; var m:TMatOfMat); // also sets all to 0
var i,j:integer;
begin
  for i:=1 to OuterNR do
    for j:=1 to OuterNC do begin
      Create1(nr,nc,m[i,j]);
      fillchar1(m[i,j],0);
    end;
end;




procedure Destroy5(OuterNR,OuterNC:integer; var m:TMatOfMat);
var i,j:integer;
begin
  for i:=1 to OuterNR do
    for j:=1 to OuterNC do
      destroy1(m[i,j]);
end;




procedure NormalizeV(var x:vector);
{x will have unit norm}
var i:integer; r2:extended;
begin
  r2:=0;
  for i:=1 to x.n do
    r2:=r2+sqr(x.v[i]);
  if r2<eps2 then
    r2:=1;
  r2:=sqrt(r2);
  for i:=1 to x.n do
    x.v[i]:=x.v[i]/r2;
end;




function NormalizeVf(var x:vector):extended; overload;
{x will have unit norm, returns norm}
var i:integer; r2:extended;
begin
  r2:=0;
  for i:=1 to x.n do
    r2:=r2+sqr(x.v[i]);
  if r2<eps2 then
    r2:=1;
  r2:=sqrt(r2);
  for i:=1 to x.n do
    x.v[i]:=x.v[i]/r2;
  NormalizeVf:=r2;
end;

procedure NormalizeV2(var x:vector);
{x will have average sum of squares = 1}
var i:integer; r2:extended;
begin
  r2:=0;
  for i:=1 to x.n do
    r2:=r2+sqr(x.v[i]);
  r2:=r2/x.n;
  if r2<eps2 then
    r2:=1;
  r2:=sqrt(r2);
  for i:=1 to x.n do
    x.v[i]:=x.v[i]/r2;
end;




procedure SubtractMean(var x:vector); overload;
{x will have zero mean}
var j:integer; r1:extended;
begin
  r1:=0;
  for j:=1 to x.n do
    r1:=r1+x.v[j];
  r1:=r1/x.n;
  for j:=1 to x.n do
    x.v[j]:=x.v[j]-r1;
end;




function Trace(var m:matrix):real;
var i:integer; r:extended;
begin
  r:=0;
  for i:=1 to m.nr do
    r:=r+m.m[i,i];
  trace:=r
end;




function UnitTrace(var m:matrix):real;
// returns tr and makes m=m/tr (forces matrix to unit trace)
var i,j:integer; r:extended;
begin
  r:=0;
  for i:=1 to m.nr do
    r:=r+m.m[i,i];
  UnitTrace:=r;
  for i:=1 to m.nr do
    for j:=1 to m.nr do
      m.m[i,j]:=m.m[i,j]/r;
end;




begin
  case sizeof(real) of
    4:begin eps1:=1.0e-7; eps2:=1.0e-38; end;
    8:begin eps1:=1.0e-15; eps2:=1.0e-308; end;
    10:begin eps1:=1.0e-19; eps2:=1.0e-4932; end;
    else
    halt;
  end;
end.

