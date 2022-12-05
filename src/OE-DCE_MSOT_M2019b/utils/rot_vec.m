function [ x1,y1,z1 ] = rot_vec( x,y,z,phi,theta )

%rotation around z
xt=cos(-phi).*x-sin(-phi).*y;
yt=sin(-phi).*x+cos(-phi).*y;
zt=z;

%rotation around y
x1=cos(theta-pi/2).*xt+sin(theta-pi/2).*zt;
y1=yt;
z1=-sin(theta-pi/2).*xt+cos(theta-pi/2).*zt;

end

