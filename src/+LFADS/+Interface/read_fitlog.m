function log = read_fitlog(filename)
% function log = lfadsi_read_fitlog(filename)
    fh = fopen(filename ,'r');
    lc=linecount(fh);
    fclose(fh);

    log = textread(filename,'%s','whitespace',',');
    log = reshape(log,[],lc)';

end

function n = linecount(fid)
n = 0;
tline = fgetl(fid);
while ischar(tline)
    tline = fgetl(fid);
    n = n+1;
end
end