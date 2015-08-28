function [A, M, B, Lambda, R, Tt, Tu, Ts, Ws, Us, Vs, EXs, EYs] = ...
    canMBPLS_compMap(Xs, Ys, varExp, p, groupLabel)
% function [A, M, B, Lambda, R, Ts, Ws, Us, Vs, EXs, EYs] = canMBPLS_compMap(Xs, Ys)
% Inputs: 
% Xs (1x number of data blocks) cell array where element s is the observed
%   data of subject 's', which has the dimension of t (time) x v (voxel)
% Ys (1x number of data blocks) cell array where element s is the observed
%   data of subject 's', which has the dimension of t (time) x z (prior 
%   information)
%
% Outputs:
% A (number of independent voxel patterns x number of voxels) Group-level
%   independent voxel patterns
% M (number of voxel patterns x number of voxel patterns) Group-level
%   mixing matrix
% B (number of group comp x number of voxels) Group-level
%   principal components given by CCA
% Lambda ( (num of subj * num subj comp) x number of group comp) CCA
%   subject loadings
% R ( (num of subj * num subj comp) x number of voxels) Subject-variability
%   noise
% Tt (number of time points x number of subj comp) mbPLS X-block superscores
% Tu (number of time points x number subj comp) mbPLS Y-block superscores
% Ts (1 x number of data blocks) cell array where element s is the
%   subject-level X-block time courses of subject 's', which has the dimension of
%   (number of time points x number of subj comp)
% Ws (1 x number of data blocks) cell array where element s is the
%   subject-level X-block spatial patterns of subject 's', which has the dimension 
%   of (number of subj comp x number of voxels)
% Us (1 x number of data blocks) cell array where element s is the
%   subject-level Y-block time courses of subject 's', which has the dimension of
%   (number of time points x number of subj comp)
% Vs (1 x number of data blocks) cell array where element s is the
%   subject-level Y-block spatial patterns of subject 's', which has the dimension 
%   of (number of subj comp x prior information)
% EXs (1 x number of data blocks) cell array where element s is the
%   subject-level X block observation noise of subject 's', which has the dimension 
%   of (number of time points x number of voxels)

if ~iscell(Xs) || ~iscell(Ys)
    error('Xs and Ys must be cell arrays');
end
if nargin < 4
    p = 0.05;
end
if nargin < 5
    groupLabel = ones(1, length(Xs));
end

numSubj = length(Xs);
[T, numVoxel] = size(Xs{1});
numPrior = size(Ys{1}, 2);
tol = 1e-12;

Ts = cell(1, numSubj);
Ws = cell(1, numSubj);
Us = cell(1, numSubj);
Vs = cell(1, numSubj);
EXs = cell(1, numSubj);
EYs = cell(1, numSubj);

% decompose each subject's data mbPLS
[Tb,Pb,Wb, Wb_reproj, Wt,Tt,Ub,Qb,Wu,Tu, ssx, ssy] = ...
    MBbiPLS2(Xs,Ys,min(15, numVoxel),tol);
% determine the optimal number of PLS factors to estimate for each subject
n_sbj = min(find(ssx{1}(:, 1) > varExp));
if isempty(n_sbj)
    warning('canMBPLS_compMap: insufficient number of PLS components. Reset to 10');
    n_sbj = 10;
end

for s = 1:numSubj
    Ts{s} = Tb{s}(:, 1:n_sbj);
    Ws{s} = Wb{s}(:, 1:n_sbj).';
    Us{s} = Ub{s}(:, 1:n_sbj);
    Vs{s} = Qb{s}(:, 1:n_sbj).';
    EXs{s} = Xs{s} - Ts{s}*Ws{s};
    EYs{s} = Ys{s} - Us{s}*Vs{s};
end

% determine the optimal number of canonical variables to keep
% keeping the first canonical variable from each PLS component
n_grp = 1;

% % % % determine resampling block size
% % % [Warfit, Aarfit, Carfit, sbc] = ...
% % %     arfit(Xs{randi(numSubj)}(:, randi(numVoxel)), 1, 20); 
% % % [junk, arOrder] = min(sbc);
% % % arOrder = 4;    % determined apriori
% % % blockSize = 2*arOrder + 1;
% % % blockIndex = (blockSize:blockSize:size(Xs{1}, 1))/blockSize;
% % % dataIndex = reshape((1:length(blockIndex)*blockSize)', blockSize, :).';
% % % % generate surrogate data
% % % % permute the time segements in X data blocks
% % % numSurr = 100;
% % % zDistr = NaN(1,numSurr*n_sbj);
% % % 
% % % for k=1:numSurr
% % %     z = randperm(length(blockIndex));
% % %     permIndex = reshape(dataIndex(z, :).', [], 1);
% % %     Xs_surr = cell(1, numSubj);
% % %     for s = 1:numSubj
% % %         Xs_surr{s} = Xs{s}(permIdex, :);
% % %     end
% % %     [Tb,Pb,Wb, Wb_reproj, Wt,Tt,Ub,Qb,Wu,Tu, ssx, ssy] = ...
% % %         MBbiPLS2(Xs_surr,Ys,min(15, numVoxel),tol);
% % %     C = corrcoef([Tt, Tu]);
% % %     numPLS = size(Tt, 2);
% % %     zDistr = [zDistr, abs(reshape(diag(C(1:numPLS, numPLS:end)), 1, []))];
% % % end
% % % zDistr = sort(zDistr, 'descend');
% % % zth = zDistr(floor(numSurr*p));
    
    


groups = unique(groupLabel);
numGroups = length(groups);
B = cell(1, n_sbj);
Lambda = cell(1, n_sbj);
R = cell(1, n_sbj);
for f = 1:n_sbj
    if numGroups == 1
        P = [];
        for s = 1:numSubj
            P = [P; Ws{s}(f, :)]; % (numSubj*n_sbj) x numVoxel
        end

        % perform generalized SVD on P
        [Upsilon,Zeta,Theta] = svd(P, 'econ');
        B{f} = Theta(:, 1:n_grp).';
        Lambda{f} = Upsilon(:, 1:n_grp) * Zeta(1:n_grp, 1:n_grp);
        R{f} = P - Lambda{f} * B{f};
    else
        B{f} = cell(1, numGroups);
        Lambda{f} = cell(1, numGroups);
        R{f} = cell(1, numGroups);
        for g = 1:numGroups
            % Apply generalized canonical correlation analysis (CCA) to estimate the 
            % group-level patterns
            % Concatenate subject-level spatial patterns Ps to form P
            P = [];
            for s = find(groupLabel == groups(g))
                P = [P; Ws{s}(f, :)]; % (numSubj*n_sbj) x numVoxel
            end

            % perform generalized SVD on P
            [Upsilon,Zeta,Theta] = svd(P, 'econ');
            B{f}{g} = Theta(:, 1:n_grp).';
            Lambda{f}{g} = Upsilon(:, 1:n_grp) * Zeta(1:n_grp, 1:n_grp);
            R{f}{g} = P - Lambda{f}{g} * B{f}{g};
        end
    end
end

if numGroups == 1
    % perform ICA on B
    % [A, M, W] = fastica(cell2mat(B.'));
     [weights,sphere,A]  = mexica(cell2mat(B.'));
     M = pinv(weights*sphere);
else
    A = cell(1, numGroups);
    M = cell(1, numGroups);
    for g = 1:numGroups
        Btemp = [];
        for f = 1:n_sbj
            Btemp = [Btemp; B{f}{g}];
        end
        [weights, sphere, A{g}] = mexica(Btemp);
        M{g} = pinv(weights * sphere);
    end
end
