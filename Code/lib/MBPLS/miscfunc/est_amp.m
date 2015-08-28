function y = est_amp( x, n, w )
% EST_AMP
%   Estimate sEMG amplititude of a channel
%
% y = est_amp( x, n, w )
%
% Input:
%   x: the sEMG signal, a vector
%   n: the order of the estimator, a positive number
%   w: the window width of the estimator, an positive integer not greater
%   than the length of x.
%
% Output:
%   y: the amplitude the sEMG signal
%
% Algorithm:
%   Three steps:
%       1) demodulate:  y = abs( x ) .^ n
%       2) smooth:      y was smoothed with a moving average window
%       3) relinearize: y = abs( x ) .^ 1/n

% demodulate
y = abs( x ) .^ n;

% smooth
y = smooth( y, w );

% relinearize
y = abs( y ) .^ ( 1 / n );

function y = smooth( x, w )
% first pass
y = ma( x, w );

% The codes below are commented because they do not work well.
% It seems that the estimated optimal window size is too large
%
% % Optimal window size, according to Clancy, IEEE Trans. on Biomed. Eng.,
% % vol 46, page 717, 1999
%
% s1 = mean( y .* y );
% ddy = diff( diff( y ) ); s2 = mean( ddy .* ddy );
%
% f = 600;
% w = ceil( f * ( ( 72 / g ) * ( s1 / s2 ) ) .^ 0.2 );
% 
% % second pass
% y = ma( x, w );

function y = ma( x, w ) % subfunction, moving average
% To get rid the boundary effect, x is expanded
% ( x0 x1 x2 .... xn-2 xn-1 xn ) =>
% ( ... x2 x1 x0 x1 x2 ... xn-2 xn-1 xn xn-1 xn-2 ... )
% After the expansion, moving average is caculated from x1 to xn.

y = zeros( size( x ) );

if mod( w, 2 )
    w1 = ( w -1 ) / 2;
    w2 = w1;
else
    w1 = w / 2;
    w2 = w1 - 1;
end

if 1 == size( x, 1 )    % row vector
    x = x( : );
end
x = [ x( ( w1 + 1 ) : -1 : 2 ); x; x( ( end - 1 ) : -1 : ( end - w2 ) ) ];

for i = 1 : length( y )
    y( i ) = mean( x( i : ( i + w - 1 ) ) );
end

