function Ri = calc_Ri(kk,N_x,N_y,thres,depth_im,wl_name)
% Ri{kk} = calc_Ri(kk,N_x,N_y,thres,depth_im,wl_name)
% ---------------------------------------------------
% Calculate the inverse wavelet transform matrix (Ri in eq. 24)
%
% - kk:         Number of wavelet base
% - N_x, N_y:   Image dimension
% - thres:      Threshold to ignore weak crosstalk (e.g. 0.2)
% - depth_im:   Number of decomposition levels (e.g. 2)
% - wl_name:    Name of the wavelet base (e.g. 'db6')

    % Lookup Table dependent on image size (needs optimisation)
    [temp,L4]=wavedec2(ones(N_x,N_y),depth_im,wl_name);
    clear temp;

    % get base and inverse for selected base
    base=return_base_arb_lev_even_db(kk,depth_im,wl_name );
    Ri_full=calc_Ri_part_arb_lev( base,N_x,N_y,L4,depth_im );

    % find important places in inverse to help ignore weak crosstalk
    places=1:size(Ri_full,2);
    energy=sum(abs(Ri_full).^2);
    max_energy=max(energy);
    places=places((energy>thres*max_energy));

    Ri = sparse(Ri_full(:,places));
end


















% calculate wavelet base function
function [base ] = return_base_arb_lev_even_db(lev_n,depth,wl_name)  %calculates the wavelet reconstruction matrix 

    [LO_D,HI_D,LO_R,HI_R] =wfilters(wl_name);
    base_cell{1}=(LO_R).'*(LO_R);
    base_cell{2}=(HI_R).'*(LO_R );
    base_cell{3}=(LO_R).'*(HI_R );
    base_cell{4}=(HI_R).'*(HI_R );

    lev_n=lev_n-1;

    lev=zeros(1,depth);
    for kk=1:depth
        lev(kk)=floor(lev_n/4^(depth-kk));
        lev_n=lev_n-4^(depth-kk)*lev(kk);
    end
    lev=lev+1;

    %%%%%%%
    % How can we calculate only one element in B every time?
    % 
    % We need to have the right building blocks. 
    % this is done by convolving the bases + down/up sampling taken into account. 
    % everything is done with jumps of 2

    base=base_cell{lev(end)};
    for ii=depth-1:-1:1
      base_spread=zeros(size(base,1)*2-1,size(base,2)*2-1);
      for kk=1:size(base,1)
         for uu=1:size(base,2)
            base_spread(kk*2-1,uu*2-1)= base(kk,uu);
         end
      end
      base=conv2(base_spread,base_cell{lev(ii)}); 
    end
end   % function

  
  
  
  
  
  
  
  
  
  
function [ Ri_full ] = calc_Ri_part_arb_lev( base,N_x,N_y,L5,depth )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


    n_xA=1:size(base,2);
    n_xA=n_xA-(length(n_xA)-2^depth );
    n_yA=1:size(base,1);
    n_yA=n_yA-(length(n_yA)-2^depth );

    u=zeros(1,L5(1)*L5(2)*N_x*N_y);

    for ii=1:length(n_xA);%length(n_xA)
        for jj=1:length(n_yA);%length(n_yA)
            start_x=n_xA(jj);
            start_y=n_yA(ii);

            base_pix=base(jj,ii);

            x_places=start_x:2^depth:start_x+(L5(1,2)-1)*2^depth;
            y_places=start_y:2^depth:N_y+(L5(1,1)-1)*2^depth;


            NN=N_x*N_y;

            [X_places,Y_places]=meshgrid(x_places,y_places);
             X_places=X_places.';
             Y_places=Y_places.';


            X_places_vec=X_places(:);
            Y_places_vec=Y_places(:);
            is_pos= and(and((X_places_vec>0),(Y_places_vec>0)),and((X_places_vec<=N_x),(Y_places_vec<=N_y)));
            %is_pos=(X_places_vec>0).*(Y_places_vec>0).*(X_places_vec<=N_x).*(Y_places_vec<=N_y);

            places=X_places_vec+N_y*(Y_places_vec-1);
            n=[1:length(places)]';
            places2=places+(n-1)*NN;
            %is_pos2=is_pos+(n-1)*NN;


            places2_pos=places2(is_pos);

            u(places2_pos)=base_pix;

        end
    end
    Ri_full=reshape(u,N_x*N_y,L5(1)*L5(2));
end      % function





