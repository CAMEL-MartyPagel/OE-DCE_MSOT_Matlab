tserial = '9967D101';

r_sensor = 0.04045;
imp_resp = [];

format long
numel = 256;
coverage = 270;
elementstep = coverage/360*2*pi/numel
firstelement = (180-coverage)/2/360*2*pi + (coverage/360*2*pi/numel/2)
angle_sensor = firstelement:elementstep:firstelement+(numel-1)*elementstep;

center =(firstelement+(angle_sensor(end)-angle_sensor(1))/2)/2/pi*360

global TFILEDIR;
save(['tdata_' tserial],'angle_sensor','r_sensor','imp_resp');

