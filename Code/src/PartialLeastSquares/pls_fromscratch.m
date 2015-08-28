function [P, Q, T, U] = pls_fromscratch(X, Y, ncomps)
%PLS_FROMSCRATCH Returns X and Y loadings from partial least squares alg.
%   The point of this function is for me to implement the basic partial
%   least squares function outlined in Geladi, Kowalski (1986) and build
%   from this to write multi-block PLS.
[n, m] = size(X);
th = 1e-10;
[ny, p] = size(Y);
if (n ~= ny)
    fprintf('Error: samples of X do not match samples of Y.\n');
    P = 0; Q = 0; return;
end

% The data must be mean-centered and variance-scaled
% X
for i=1:m
    % subtract mean
    X(:,i) = X(:,i) - mean(X(:,i))*ones(n,1);
    % variance scaled
    X(:,i) = X(:,i) / std(X(:,i));
end
% Y
for i=1:p
    % subtract mean
    Y(:,i) = Y(:,i) - mean(Y(:,i))*ones(n,1);
    % variance scaled
    Y(:,i) = Y(:,i) / std(Y(:,i));
end

P = zeros(m,ncomps);
Q = zeros(p,ncomps);
T = zeros(n,ncomps);
U = zeros(n,ncomps);
E = X;
F = Y;

for h=1:ncomps
    X = E;
    Y = F;
    j = randi(p);
    u = Y(:,j);                  % (1)
    t_old = -10*ones(n, 1);
    t = zeros(n,1);
    ind = 1;
    while (norm(t - t_old) > th) % (8)
        fprintf('Iteration %d\n', ind);
        % In the X block
        w = (u'*X / (u'*u))';    % (2)
        w = w / norm(w);         % (3)
        t_old = t;
        t = X*w / (w'*w);        % (4)
        % In the Y block
        q = (t'*Y / (t'*t))';    % (5)
        q = q / norm(q);         % (6)
        u = Y*q / (q'*q);        % (7)
        ind = ind + 1;
    end
    % Calculate X loadings and rescale the scores and weights accordingly
    ph = (t'*X / (t'*t))';       % (9)  (the variable ph is used to distingish from p number of output vars)
    C = norm(ph);
    ph = ph / C;                 % (10)
    t = t / C;                   % (11)
    w = w / C;                   % (12)
    P(:,h) = ph;
    Q(:,h) = q;
    T(:,h) = t;
    U(:,h) = u;
    
    % Find the regression coefficient b for the inner relation
    b = u'*t / (t'*t);           % (13)
    
    % Calculate the residuals
    E = E - t*ph';               % (14)
    F = F - b*t*q';              % (15)
end
    

end

