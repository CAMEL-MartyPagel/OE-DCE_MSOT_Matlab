function [zpositions,timestamps, wavelengths] = get_datainfo(datainfo)
%GET_DATAINFO Extract Information
%  [zpos, timestamps] = get_datainfo(datainfo) 
%     Extracts information about contained z-slices and timepoints of
%     measurements from the datainfo struct that is returned by loadData.
%
%  timestamps contains the timestamp of the first wavelength in seconds
%
%  See also LOADDATA


wavelengths = datainfo.measurement_desc.wavelengths(:,2)*1e9;

frames = datainfo.measurement_desc.projections.msot_frame;
zpositions = zeros(size(frames));
for i = 1:size(frames,1)
    zpositions(i) = frames{i,1}.zpos;
end;
clear i ;
zpositions = unique(zpositions)*1e3;

timestamps = zeros(uint8(numel(frames)/numel(wavelengths)/numel(zpositions)),1);
for i = 1:size(frames,1)/numel(wavelengths)/numel(zpositions)
    timestamps(i) = frames{(i-1)*numel(wavelengths)*numel(zpositions)+1,1}.timestamp*1e-3;
end;
