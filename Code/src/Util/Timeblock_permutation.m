function [Y] = Timeblock_permutation(X, blocksize, permutation)
%TIMEBLOCK_PERMUTATION Summary of this function goes here
%   Detailed explanation goes here
    nblocks = ceil(size(X,1)/blocksize);
    if (nblocks ~= length(permutation))
        fprintf('Error: permutation length %d does not match number of block %d.\n', ...
                length(permutation), nblocks);
    end
    
    % Store blocks
    blocks = cell(nblocks, 1);
    for i=1:(nblocks-1)
        blocks{i}.data = X((i-1)*blocksize+1:i*blocksize, :);
        blocks{i}.len = blocksize;
    end
    blocks{nblocks}.data = X((nblocks-1)*blocksize+1:end,:);
    blocks{nblocks}.len = size(blocks{nblocks}.data, 1);
    
    % Permute blocks into Y
    Y = zeros(size(X));
    blockstart = 1;
    for i=1:nblocks
        j = permutation(i);
        block_data = blocks{j}.data;
        block_len = blocks{j}.len;
        Y(blockstart:(blockstart+block_len-1),:) = block_data;
        blockstart = blockstart + block_len;
    end
    

end

