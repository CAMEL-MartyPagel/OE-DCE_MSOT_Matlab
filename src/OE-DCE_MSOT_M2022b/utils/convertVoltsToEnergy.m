function en = convertVoltsToEnergy(corr,wl,val)
% CONVERTVOLTSTOENERGY
% res = CONVERTVOLTSTOENERGY(correctionTable,wavelength,diodeReadout)
%
% Converts diode readout from datainfo.ScanFrame(x).DiodeReadout to Laser
% Energy in mJ using the correction table supplied.
%
% See also: readVoltsToEnergy


wlind = find(corr(:,1) == wl);      % find wavelength index
pol = corr(wlind,2:end);            % find coefficients of polynomial

% calculate energy by evaluating 5th order polynomial
en = polyval(fliplr(pol),val);