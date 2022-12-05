function [ A_out ] = transform_t2wl( A, n_t,n_proj,depth,wl_name ) %Performs time-proj 2d-wlet directly on the matrix
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
temp=calc_db2_wlet_dec2_pack_arb_lev_even_db(full(reshape(A(:,1),n_t,n_proj)),depth,wl_name).';
A_out=(zeros(length(temp),size(A,2)));
A_out(:,1) = temp;

%if you want to conserve memory use A_out=sparse(zeros(length(temp),size(A,2))); 
%this option is slower though

for kk=2:size(A,2)  %this calculates the multipication of C_doulbe with A_mat, more accurate than the matrix
  A_out(:,kk)=calc_db2_wlet_dec2_pack_arb_lev_even_db(full(reshape(A(:,kk),n_t,n_proj)),depth,wl_name).';
%   if mod(kk,100)==1
%       kk
%   end
end
A_out=sparse(A_out);

end





function [ output ] = calc_db2_wlet_dec2_pack_arb_lev_even_db( im,depth,wl_name) % decomposes the 2D signalto wavepacket coefficients 
%depth must be >=1

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%len=length(w_coeff)/4;

[sig_temp,L]=wavedec2(im,1,wl_name);

if depth>1
sig_temp1=calc_db2_wlet_dec2_pack_arb_lev_even_db(reshape(sig_temp(1:end/4),L(1,1),L(1,2)) ,depth-1,wl_name);
sig_temp2=calc_db2_wlet_dec2_pack_arb_lev_even_db(reshape(sig_temp(end/4+1:2*end/4),L(1,1),L(1,2)) ,depth-1,wl_name);
sig_temp3=calc_db2_wlet_dec2_pack_arb_lev_even_db(reshape(sig_temp(2*end/4+1:3*end/4),L(1,1),L(1,2)) ,depth-1,wl_name);
sig_temp4=calc_db2_wlet_dec2_pack_arb_lev_even_db(reshape(sig_temp(3*end/4+1:4*end/4),L(1,1),L(1,2)) ,depth-1,wl_name);

    
    
%     
% [sig_temp1,LL]=wavedec2(reshape(sig_temp(1:end/4),L(1,1),L(1,2)),1,'db2');
% [sig_temp2,LL]=wavedec2(reshape(sig_temp(end/4+1:2*end/4),L(1,1),L(1,2)),1,'db2');
% [sig_temp3,LL]=wavedec2(reshape(sig_temp(2*end/4+1:3*end/4),L(1,1),L(1,2)),1,'db2');
% [sig_temp4,LL]=wavedec2(reshape(sig_temp(3*end/4+1:4*end/4),L(1,1),L(1,2)),1,'db2');
% 



output=[sig_temp1,sig_temp2,sig_temp3,sig_temp4];

else
    
output=sig_temp;

end
end

