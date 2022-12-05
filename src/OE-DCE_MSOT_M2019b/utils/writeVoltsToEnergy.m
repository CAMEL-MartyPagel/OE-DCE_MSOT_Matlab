function writeVoltsToEnergy(corr,fname)
% WRITEVOLTSTOENERGY
%   WRITEVOLTSTOENERGY(correctionTable,fileName)
%
% Writes Per Pulse correction Table to the binary file specified in
% fileName.
%
% Format is 301 wavelenghts by 7 columns (column 1 is wavelenght)
%
% See also: readVoltsToEnergy

if (size(corr) ~= [301 7])
    error('Incorrect table size - should be 301x7');
end

% open and load data
FID = fopen(fname,'w');
data = fwrite(FID,corr,'double');
fclose(FID);
