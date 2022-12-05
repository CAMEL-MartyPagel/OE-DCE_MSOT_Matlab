
function [SO2] = msp_oxymap(Rec,Hb,HbO2,T1L,method,varargin)
% function oxy_map calculates oxygen saturation and displays it in
% colorscale from red to blue through white 
% normalization and flexible thresholding does not ensure quantitative results

% INPUT: Rec - reconstructed images for several wavelengthes (1 run, 1 position, 1 frame)
%        Hb - matrix of deoxy hemoglobin signal intensities
%        HBO2 - matrix of oxy hemoglobin signal intensities
%        ThbL - low thresholding for deoxy hemoglobin
%        Thbo2L - low thresholding for oxy hemoglobin

% note: increasing both thresholds leads to highlighting regions only of very high or
% very low saturation

%close all;


% load('Rec.mat');
% load('Hb.mat');
% load('HbO2.mat');


if numel(varargin) == 1
    mask = varargin{1};
end

[cmap]=buildcmap('bwr');

Rec = squeeze(Rec);

% settings
filter = 0;
norm = 0;
imod = 0;
mask_flag = 0;
overlay = 0;


switch method 
%%
    case 'SO2'

   
       SO2 = zeros(size(Hb));
        
       % define invalid pixels where Hb OR HbO2 are less or equal zero
       pixInvalid = (Hb<=0) | (HbO2<=0);
        
       % discard negatives, compute HbT
       Hb = max(Hb,0);
       HbO2 = max(HbO2,0);
       HbT = HbO2 + Hb; 

       % compute SO2, set invalid pixels to NaN
       SO2 = HbO2./HbT; % modified EM 20161011
       SO2(pixInvalid) = NaN;
        
       % set invalid pixels to NaN
       T = T1L;
       SO2(HbT<=T) = NaN;  
    
   
% % Gaussian smoothing
% sigma = 0.5; M_g  = fspecial('gaussian',ceil(3*sqrt(sigma)),sigma); %%Gaussian Mask 
% SO2 = imfilter(SO2,M_g,'replicate');

case 'SO2_Adrian'

   Hb = max(Hb,0);
   HbO2 = max(HbO2,0);
        
   HbT = HbO2 + Hb;
   SO2= zeros(size(HbT));
   
   SO2(HbT>0) = HbO2(HbT>0)./(HbT(HbT>0));
   
   T = T1L*max(HbT(:)); % from relative to absolute thresholds
      
   SO2(HbT<T) = 0; % thresholding of the regions with negative values only

case 'difference'
 
     numWL1 = varargin{1};
     numWL2 = varargin{2};
     
     Rec1 = Rec(:,:,numWL1);
     Rec2 = Rec(:,:,numWL2);
    
     T1L = T1L*max(Rec1(:)); % from relative to absolute thresholds
     T2L = T2L*max(Rec2(:));
      

     
     Rec1(Rec1 < T1L) = T1L; % thresholding
     Rec2(Rec2 < T2L) = T2L;
     
     SO2 = (Rec2 - Rec1);  % difference
     SO2 = SO2 - min(SO2(:));
     
case 'ratio'
   
     numWL1 = varargin{1};
     numWL2 = varargin{2};
     limL = varargin{3};
     limH = varargin{4};
     
     Rec1 = Rec(:,:,numWL1);
     Rec2 = Rec(:,:,numWL2);
    
     T1L = T1L*max(Rec1(:)); % from relative to absolute thresholds
     T2L = T2L*max(Rec2(:));
        
  
     
     Rec1(Rec1 < T1L) = T1L; % thresholding
     Rec2(Rec2 < T2L) = T2L;
     
       
     SO2 = (Rec2 - Rec1)./Rec1; % ratio
     SO2 = SO2 - min(SO2(:));
     SO2(SO2>1) = 1;
        
end

    if norm == 1;
        % modulate with intensity
         if imod == 1; 
           bg = Rec(:,:,:,1);
           bg = mat2gray(bg);
           SO2(:,:,:) = SO2(:,:,:).*bg.^(1/2);
        end;
        
        SO2(:,:,:) = newrange(SO2(:,:,:),[0 1]);
    end;
  
    
    if filter == 1;
   
        SO2(:,:,:) = medfilt2(SO2(:,:,:),[3 3]);
    
    end

    
%% 


fg = SO2(:,:,:,1); % foreground

%% Display options
if mask_flag == 1 && overlay == 0
    
   figure;imagesc(fg);  colormap(cmap); axis square; axis off;  colorbar; title(sprintf(method));
   caxis([0.0 1]);
   alpha_bg = 1.0;
   
   alpha_fg = ~isnan(SO2);
   alpha_fg(isnan(SO2)) = 0;
   alpha_fg = alpha_fg.*mask;
   bg = zeros(size(SO2));
   figure;image_overlay_black(bg,fg,alpha_bg,alpha_fg,cmap);title('SO2'); % overlay image anatomy and oxy map 
 
%%    
elseif mask_flag == 0 && overlay == 1  
   bg = Rec(:,:,1); % background
   bg = mat2gray(bg);
    
   alpha_bg = 1.0;
   
   % either
   alpha_fg = ~isnan(SO2);
   alpha_fg(isnan(SO2)) = 0;
   % or  
   %alpha_fg = bg.^(1/2);
    
   %Gaussian smoothing
   sigma = 0.7; M_g  = fspecial('gaussian',ceil(6*sqrt(sigma)),sigma); %%Gaussian Mask 
   alpha_fg = imfilter(alpha_fg.^2,M_g,'replicate');
   
   image_overlay(bg,fg,alpha_bg,alpha_fg,cmap);title('SO2'); % overlay image anatomy and oxy map 

%%   
elseif mask_flag == 0 && overlay == 0 

%    %imagesc(fg);  colormap(cmap); axis square; axis off;  colorbar; title(sprintf(method));
%    %caxis([0.0 1]);
%    alpha_bg = 1.0;
%    
%    alpha_fg = ~isnan(SO2);
%    alpha_fg(SO2==NaN) = 0;
%    %alpha_fg = alpha_fg.*mask;
%    bg = zeros(size(SO2));
%    image_overlay_black(bg,fg,alpha_bg,alpha_fg,cmap); % overlay image anatomy and oxy map 

%%
elseif mask_flag == 1 && overlay == 1 
    bg = mat2gray(Rec(:,:,1)); % background 
    %Gaussian smoothing
    %  sigma = 0.8; M_g  = fspecial('gaussian',ceil(3*sqrt(sigma)),sigma); %%Gaussian Mask 
    %  mask = imfilter(mask,M_g,'replicate');


 fg = fg.*mask;
 alpha_fg = mask.*bg.^(1/2); 
 alpha_fg(isnan(fg)) = 0;
 
% alpha_fg = fg;
% alpha_fg= alpha_correct(alpha_fg).*mask;

 alpha_bg = 1.0;
   
 %Gaussian smoothing
 sigma = 0.5; M_g  = fspecial('gaussian',ceil(3*sqrt(sigma)),sigma); %%Gaussian Mask 
 %alpha_fg = imfilter(alpha_fg,M_g,'replicate');
 %fg = imfilter(fg,M_g,'replicate');
  
 image_overlay(bg,fg,alpha_bg,alpha_fg,cmap); % overlay image anatomy and oxy map 
end;    

end