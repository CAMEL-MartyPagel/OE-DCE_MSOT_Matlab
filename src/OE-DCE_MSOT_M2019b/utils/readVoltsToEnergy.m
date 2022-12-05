function data = readVoltsToEnergy(fname)

numwl = 301;            % number of wavelengths (default: 680-980 == 301)
numcol = 7;             % number of columns (wl, coeff for 5th grade)
data = [];              % return vector

if isempty(fname)
    fname = 'VoltToEnergyVsWavelength.bin';
end

% check if file exists
cfiles = dir(fname);
if ~numel(cfiles)
    error(['File ' fname ' does not exist']);
end

% check file size
if cfiles(1).bytes ~= 8*numwl*numcol
    error(sprintf('Incorrect file size (is: %i bytes, should be: %i bytes)',cfiles(1).bytes,8*numwl*numcol));
end

% open and load data
FID = fopen(fname,'r');
data = fread(FID,[numwl numcol],'double');
fclose(FID);

