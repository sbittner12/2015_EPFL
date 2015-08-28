function [chlocs_new] = convertToXYZ(chfname, channels)
%%CONVERTTOXYZ - converts .ced file to .xyz file for Cartool given channels
%   convertToXYZ(chfname, channels)
%
%   PARAMETERS:
%     chfname - string : filename of channel locations with .ced file ext
%     channels - cell (1 x nchans) of strings : electrode ids
%
%   RETURNS:
%     none
%
%   @ 2015 Sean Bittner sbittner@andrew.cmu.edu   
%
chlocs = readlocs(chfname);
nchans = length(channels);
channels_to_keep = zeros(1,nchans);
j = 1;
for i=1:64
    if (nnz(ismember(channels, chlocs(i).labels)))
        channels_to_keep(j) = i;
        j = j + 1;
    end
end

chlocs_new = chlocs(channels_to_keep);

output=fopen(sprintf('electrodes_%d.xyz', nchans),'w');
fprintf(output,'%u\t%u',nchans,1);
for nn=1:nchans
    fprintf(output,'\n%f\t%f\t%f\t%s',-chlocs_new(nn).X,...
        -chlocs_new(nn).Y, chlocs_new(nn).Z,...
        chlocs_new(nn).labels);
end
fclose(output);

end

