function [sO2map]= msp_nn_sO2(Recon,WL)
% MSP_NN_SO2 Unmixing for sO2 using a neural network
%
% Last dimension must be wavelengths, returns unmixed images in same
% dimensionality 
%
% Parameters:
% - Recon: Reconstructed images, 3-XD matrix with wavelength as last dim
% - WL: Vector with used wavelengths
%
% Return Values:
% - sO2map

% add current folder to python path
insert(py.sys.path,int32(0),'');

% import python function
import py.msp_nn_sO2_interface.msp_nn_sO2_py

% unmix
sO2_py=msp_nn_sO2_py(Recon(:)',WL);

% convert pyhton array to matlab array
sO2_mat = double(py.array.array('d',py.numpy.nditer(sO2_py))); 

% reshape
svec=size(Recon);
sO2map=reshape(sO2_mat,svec(1:(end-1)));

end
