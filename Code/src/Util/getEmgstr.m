function emgstr = getEmgstr(exptype)
%GETEMGSTR - Provides filename string for opening an EMG file
%   Based on the experiment type, the proper string is provided for
%   accessing the corresponding EMG file.
%
%  PARAMETERS
%    exptype - string for experiment type based on exoskeleton function
%      Ex: 'Healthy', 'IAF', 'LT', 'MJ', or 'OAF'
%
%  RETURNS
%    emgstr - string that is part of the corresponding EMG filename
%
%  @ 2015 Sean Bittner   sbittner@andrew.cmu.edu
%
    switch exptype
        case 'Healthy'
            emgstr = 'EE_healthy';
        case 'IAF'
            emgstr = 'IAF';
        case 'LT'
            emgstr = 'EE_linear';
        case 'MJ'
            emgstr = 'EE_MJ';
        case 'OAF'
            emgstr = 'OAF';
        otherwise
            fprintf('Error: invalid exptype.\n');
            emgstr = '';
    end

end

