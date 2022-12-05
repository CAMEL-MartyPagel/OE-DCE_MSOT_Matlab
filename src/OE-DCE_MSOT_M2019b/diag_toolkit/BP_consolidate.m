clear all;

load('01_20deg','l_fwhm','r_fwhm','l_imax','r_imax','l_center','r_center');
l_fwhm01 = l_fwhm; r_fwhm01 = r_fwhm;
l_imax01 = l_imax; r_imax01 = r_imax;
l_center01 = l_center(l_imax); r_center01 = r_center(r_imax);
load('02_24deg','l_fwhm','r_fwhm','l_imax','r_imax','l_center','r_center');
l_fwhm02 = l_fwhm; r_fwhm02 = r_fwhm;
l_imax02 = l_imax; r_imax02 = r_imax;
l_center02 = l_center(l_imax); r_center02 = r_center(r_imax);
load('03_28deg_25avg','l_fwhm','r_fwhm','l_imax','r_imax','l_center','r_center');
l_fwhm03 = l_fwhm; r_fwhm03 = r_fwhm;
l_imax03 = l_imax; r_imax03 = r_imax;
l_center03 = l_center(l_imax); r_center03 = r_center(r_imax);
load('04_33deg_25avg','l_fwhm','r_fwhm','l_imax','r_imax','l_center','r_center');
l_fwhm04 = l_fwhm; r_fwhm04 = r_fwhm;
l_imax04 = l_imax; r_imax04 = r_imax;
l_center04 = l_center(l_imax); r_center04 = r_center(r_imax);
load('05_34deg_25avg','l_fwhm','r_fwhm','l_imax','r_imax','l_center','r_center');
l_fwhm05 = l_fwhm; r_fwhm05 = r_fwhm;
l_imax05 = l_imax; r_imax05 = r_imax;
l_center05 = l_center(l_imax); r_center05 = r_center(r_imax);
load('06_34deg','l_fwhm','r_fwhm','l_imax','r_imax','l_center','r_center');
l_fwhm06 = l_fwhm; r_fwhm06 = r_fwhm;
l_imax06 = l_imax; r_imax06 = r_imax;
l_center06 = l_center(l_imax); r_center06 = r_center(r_imax);
load('07_34deg','l_fwhm','r_fwhm','l_tdist','r_tdist','l_imax','r_imax','r_sensor','l_center','r_center');
l_fwhm07 = l_fwhm; r_fwhm07 = r_fwhm;
l_imax07 = l_imax; r_imax07 = r_imax;
l_center07 = l_center(l_imax); r_center07 = r_center(r_imax);

save 'BP_fwhm_data';

%% load and save simulations

clear all;

load('D:\work\2011-11_new_transducer\simulations\MSOTII_3D_field_hires_intEIR_ele14_azi0k7_COR40mm_Rcurv40mm_256ele','field3d','x','y');
dat40s = squeeze(field3d(41,:,:));
simx = x; simy = y;
load('D:\work\2011-11_new_transducer\simulations\MSOTII_3D_field_hires_intEIR_ele14_azi0k7_COR40mm_Rcurv38mm_256ele','field3d');
dat38s = squeeze(field3d(41,:,:));
load('D:\work\2011-11_new_transducer\simulations\MSOTII_3D_field_hires_intEIR_ele14_azi0k7_COR40mm_Rcurv37mm_256ele','field3d');
dat37s = squeeze(field3d(41,:,:));
load('D:\work\2011-11_new_transducer\simulations\MSOTII_3D_field_hires_intEIR_ele14_azi0k7_COR40mm_Rcurv36mm_256ele','field3d');
dat36s = squeeze(field3d(41,:,:));
load('D:\work\2011-11_new_transducer\simulations\MSOTII_3D_field_hires_intEIR_ele14_azi0k7_COR40mm_Rcurv35mm_256ele','field3d');
dat35s = squeeze(field3d(41,:,:));
load('D:\work\2011-11_new_transducer\simulations\MSOTII_3D_field_hires_intEIR_ele14_azi0k7_COR40mm_Rcurv42mm_256ele','field3d');
dat42s = squeeze(field3d(41,:,:));
clear field3d x;

for i = 1:size(dat40s,1)
    line = dat40s(i,:);
    fwhm40(i) = nnz(line > max(line(:))*0.5);
    line = dat35s(i,:);
    fwhm35(i) = nnz(line > max(line(:))*0.5);
    line = dat36s(i,:);
    fwhm36(i) = nnz(line > max(line(:))*0.5);
    line = dat37s(i,:);
    fwhm37(i) = nnz(line > max(line(:))*0.5);
    line = dat38s(i,:);
    fwhm38(i) = nnz(line > max(line(:))*0.5);
    line = dat42s(i,:);
    fwhm42(i) = nnz(line > max(line(:))*0.5);
end
clear i y line;
clear dat35s dat36s dat37s dat38s dat40s dat42s;
save 'simulations';


