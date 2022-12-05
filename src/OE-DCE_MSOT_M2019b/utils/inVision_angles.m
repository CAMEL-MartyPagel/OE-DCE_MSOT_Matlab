format long
coverage = 270;
nel = 512;

elementstep = coverage/360*2*pi/nel
firstelement = (180-coverage)/2/360*2*pi + (coverage/360*2*pi/nel/2);
angle_sensor = firstelement:elementstep:firstelement+(nel-1)*elementstep;
[angle_sensor(1) angle_sensor(end)]
center =(firstelement+(angle_sensor(end)-angle_sensor(1))/2)/2/pi*360