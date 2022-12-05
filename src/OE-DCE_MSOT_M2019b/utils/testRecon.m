n = 100;
filter_f = [50e3 7e6];
proj = 256+255;
iter = 0;
timeres = 3;
c = 1530;
limits = [r_sensor-0.02 r_sensor+0.02];
fs = 4e7;
r_sensor = 0.0405;
load('C:\Code\tdata_256A101.mat' ,'angle_sensor');

BP = reconWrapper([],sig,[],n,proj,r_sensor,angle_sensor,c,'direct',...
    fs,limits,timeres,filter_f,iter,1,4,1);