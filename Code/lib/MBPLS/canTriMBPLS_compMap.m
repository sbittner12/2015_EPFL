function [A, M, B, Lambda, R, Tt, Tu, W1b, Ts, Ws, Us, Vs, EXs, EYs, sigComp,ssx, ssy] = ...
    canTriMBPLS_compMap(Xs, Ys, varExp, p, groupLabel)
% function [A, M, B, Lambda, R, Ts, Ws, Us, Vs, EXs, EYs] = canTriMBPLS_compMap(Xs, Ys)
% Inputs: 
% Xs (1x number of data blocks) cell array where element s is the observed
%   data of subject 's', which has the dimension of t (time) x 
%   f (frequency) x v (voxel)
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
% W1b 
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

sigComp = [];
numSubj = length(Xs);
[T, numVoxel] = size(Xs{1});
numPrior = size(Ys{1}, 2);
tol = 1e-12;
defaultNumPLS = 10;

Ts = cell(1, numSubj);
Ws = cell(1, numSubj);
Us = cell(1, numSubj);
Vs = cell(1, numSubj);
EXs = cell(1, numSubj);
EYs = cell(1, numSubj);

% decompose each subject's data mbPLS
[Tb,W1b,W2b,Wt,Tt,Ub,Qb,Wu,Tu, ssx, ssy] = ...
    MBtriPLS2_cmmnW1b(Xs,Ys,min(defaultNumPLS, numVoxel),tol);
% determine the optimal number of PLS factors to estimate for each subject
n_sbj = min(find(ssx{1}(:, 1) > varExp));
if isempty(n_sbj)
    warning('canTriMBPLS_compMap: insufficient number of PLS components. Reset to 10');
    n_sbj = defaultNumPLS;
end

W1b = W1b{1}(:, 1:n_sbj);
for s = 1:numSubj
    Ts{s} = Tb{s}(:, 1:n_sbj);
    Ws{s} = W2b{s}(:, 1:n_sbj).';
    Us{s} = Ub{s}(:, 1:n_sbj);
    Vs{s} = Qb{s}(:, 1:n_sbj).';
%     EXs{s} = Xs{s} - reshape(Tt(:, 1:n_sbj) * kron(W2b{s}(:, 1:n_sbj), W1b(:, 1:n_sbj)).', size(Xs{s}));
%     EYs{s} = Ys{s} - Us{s}*Vs{s};
end

% determine the optimal number of canonical variables to keep
n_grp = 1;

% % % % determine resampling block size
% % % [Warfit, Aarfit, Carfit, sbc] = ...
% % %     arfit(Xs{randi(numSubj)}(:, randi(numVoxel)), 1, 20); 
% % % [junk, arOrder] = min(sbc);
% % % arOrder = 7;    % determined apriori
% % % blockSize = 2*arOrder + 1;
% % % blockIndex = (blockSize:blockSize:size(Xs{1}, 1))/blockSize;
% % % dataIndex = reshape((1:length(blockIndex)*blockSize)', blockSize, []).';
% % % % generate surrogate data
% % % % permute the time segements in X data blocks
% % % numSurr = 10;
% % % zDistr = [];
% % % 
% % % for k=1:numSurr
% % %     z = randperm(length(blockIndex));
% % %     permIndex = reshape(dataIndex(z, :).', [], 1);
% % %     Xs_surr = cell(1, numSubj);
% % %     Ys_surr = cell(1, numSubj);
% % %     for s = 1:numSubj
% % %         Xs_surr{s} = Xs{s}(permIndex, :, :);
% % %         Ys_surr{s} = Ys{s}(1:length(blockIndex)*blockSize, :);
% % %     end
% % %     [Tb_surr,W1b_surr,W2b_surr,Wt_surr,Tt_surr,Ub_surr,Qb_surr,...
% % %         Wu_surr,Tu_surr] = ...
% % %         MBtriPLS2_cmmnW1b(Xs_surr,Ys_surr,n_sbj,tol);
% % %     C = corrcoef([Tt_surr, Tu_surr]);
% % %     numPLS = size(Tt_surr, 2);
% % %     zDistr = [zDistr, abs(reshape(diag(C(1:numPLS, numPLS+1:end)), 1, []))];
% % % end
% % % zDistr = sort(zDistr, 'descend');
% % % zth = zDistr(ceil(length(zDistr)*p));
    
numPLS = n_sbj;
C = corrcoef([Tt, Tu]);
cc = abs(diag(C(1:numPLS, numPLS+1:end)));
sigComp = find(cc >= 0);

groups = unique(groupLabel);
numGroups = length(groups);
B = cell(1, n_sbj);
Lambda = cell(1, n_sbj);
R = cell(1, n_sbj);
for f = 1:n_sbj
    % Apply generalized canonical correlation analysis (CCA) to estimate the 
    % group-level patterns
    % Concatenate subject-level spatial patterns Ps to form P
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
            
% % %             y = cell(1, sum(groupLabel == groups(g)) );
% % %             counter=1;
% % %             for s = find(groupLabel == groups(g))
% % %                 y{counter} = Ws{s}(f, :);
% % %                 counter = counter + 1;
% % %             end
% % %             [Btemp] = mcca_ssqcor(y, n_grp);
                
            % perform generalized SVD on P
            [Upsilon,Zeta,Theta] = svd(P, 'econ');
            B{f}{g} = Theta(:, 1:n_grp).';
            Lambda{f}{g} = Upsilon(:, 1:n_grp) * Zeta(1:n_grp, 1:n_grp);
            R{f}{g} = P - Lambda{f}{g} * B{f}{g};
        end
    end
end

% perform ICA on B
% [A, M, W] = fastica(cell2mat(B.'));
if numGroups ==  1
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
        try
            [weights, sphere, A{g}] = mexica(Btemp);
            M{g} = pinv(weights * sphere);
        catch
            A{g} = [];
            M{g} = [];
        end
    end
end
    
