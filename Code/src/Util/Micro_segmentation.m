Microstate extraction

data = averaged epochs;

pop_writeeeg(data,['subject_' subject '_reaching.bdf'],'TYPE','BDF');
output=fopen(['subject_' subject '_reaching_electrodes.xyz'],'w');
len=length(data.chanlocs);
fprintf(output,'%u\t%u',len,1);
for nn=1:len
    fprintf(output,'\n%f\t%f\t%f\t%s',-data.chanlocs(nn).X,...
        -data.chanlocs(nn).Y,data.chanlocs(nn).Z,...
        data.chanlocs(nn).labels);
end
fclose(output);

Segmentation

micro = load('.ep');

data = concatenated epochs;

for m = 1:size(micro)
    for i = 1:time
        corr(m,i) = corr2(data(i),micro(m,:));
    end
end


