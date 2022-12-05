cfiles = rdir('**\Scan_*.bin');
swpos1 = 129;
swpos2 = 161;
swcount = 32;
nchannels = 256;
chbytes = 2030*2;   % 2 bytes (16bit) per sample, 2030 samples

for j = 1:numel(cfiles)
    fn = cfiles(j).name;
    fsize = cfiles(j).bytes;
    fprintf('File: %s\n',fn);

    % Backup
    bakfn = strrep(fn,'.bin','.bak');
    if (~exist(bakfn,'file'))
        fprintf('  Creating backup copy of binary data file\n');
        copyfile(fn,bakfn);
    end;        
        
    % check if already converted
    donefn = strrep(fn,'.bin','.done');
    if exist(donefn,'file')
        warning('  A .done file exists in this folder. Skipping since channels have been switched already\n');
        continue;
    else 
        DONEID = fopen(donefn,'w');
        fwrite(DONEID,'1');
        fclose(DONEID);
    end

    % actually convert
    fprintf('  Writing new data file\n');
    FID = fopen(bakfn,'r');
    NID = fopen(fn,'w');
    while (1)
        % read
        tmp = fread(FID,[2030 nchannels],'uint16');
        if (feof(FID)) break; end;

        fwrite(NID,tmp(:,[1:swpos1-1 swpos2:swpos2+swcount-1 swpos1:swpos1+swcount-1 swpos2+swcount:nchannels]),'uint16');
    end
    fclose(FID);
    fclose(NID);
    fprintf('OK\n');
    
end