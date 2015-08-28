function eventstr = getEventstr(exptype)
%GETEVENTSTR - Provides filename string for opening an event file
%   Based on the experiment type, the proper string is provided for
%   accessing the corresponding event file.
%
%  PARAMETERS
%    exptype - string for experiment type based on exoskeleton function
%      Ex: 'Healthy', 'IAF', 'LT', 'MJ', or 'OAF'
%
%  RETURNS
%    eventstr - string that is part of the corresponding event filename
%
%  @ 2015 Sean Bittner   sbittner@andrew.cmu.edu
%

    switch exptype
        case 'Healthy'
            eventstr = 'HM';
        case 'IAF'
            eventstr = 'IAF';
        case 'LT'
            eventstr = 'LT';
        case 'MJ'
            eventstr = 'MJ';
        case 'OAF'
            eventstr = 'OAF';
        otherwise
            fprintf('Error: invalid exptype.\n');
            eventstr = '';
    end

end

