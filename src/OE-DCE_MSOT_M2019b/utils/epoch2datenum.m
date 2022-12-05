function ts = epoch2datenum(epoch)
% correct with correction for local time zone
corr = 1;
% in case we are dealing with epoch in milliseconds
if epoch > 2e9, corr = 1000; end
ts = (epoch/8.64e4/corr+datenum(1970,1,1,0,0,0)-(java.util.Date().getTimezoneOffset/60)/24);