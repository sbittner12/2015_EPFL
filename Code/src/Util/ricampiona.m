function output = ricampiona(input,n)

c = length(input);
x=1:c;
xx=0:(c/n):c;
output = spline(x, input, xx);

output(1)=[];

%Error : [X,Y,SIZEY] = CHCKXY(X,Y) checks the data sites X and corresponding data
%   values Y, making certain that there are exactly as many sites as values,
%   that no two data sites are the same, removing any data points that involve 
%   NaNs, reordering the sites if necessary to ensure that X is a strictly
%   increasing row vector and reordering the data values correspondingly,
%   and reshaping Y if necessary to make sure that it is a matrix, with Y(:,j)
%   the data value corresponding to the data site X(j), and with SIZEY the
%   actual dimensions of the given values. 
%   This call to CHCKXY is suitable for PCHIP.
%
%   [X,Y,SIZEY,ENDSLOPES] = CHCKXY(X,Y) also considers the possibility that
%   there are two more data values than there are data sites.
%   If there are, then the first and the last data value are removed from Y
%   and returned separately as ENDSLOPES. Otherwise, an empty ENDSLOPES is
%   returned.  This call to CHCKXY is suitable for SPLINE.

% x = 0:10;
% y = sin(x);
% xx = 0:.25:10;
% yy = spline(x,y,xx);
% plot(x,y,'o',xx,yy)