function [] = run_DCE_MSOT_MSP (msot_file, saving_DIR, segm)

% Add MSOT MAtlab library to path 
% This library was created by iThera
addpath( genpath('com.itheramedical.msotlib_beta_rev723') )
javaaddpath com.itheramedical.msotlib_beta_rev723\MSOTBeans\xbean.jar %%AK
javaaddpath com.itheramedical.msotlib_beta_rev723\MSOTBeans\msotbeans.jar %%AK

%%
%   loadMSOTMsp  Load MSOT Specrally Unmixed Data
%   [R wls zpos ts datainfo] = LOADMSOTMSP()       
%           opens Open File Dialog to choose file, then loads the first
%           Msp
%   [R wls zpos ts datainfo] = LOADMSOTMSP(file)  
%           opens the specified file, then loads the first msp
%   [R wls zpos ts datainfo] = LOADMSOTMSP(datainfo)   
%           opens the specified file, then loads the first msp
%   [R wls zpos ts datainfo] = LOADMSOTMSP(datainfo,mspNum)   
%           opens the specified file, then loads the specified msp
%   [R wls zpos ts datainfo] = LOADMSOTMSP(datainfo,mspNum,selMat)   
%           opens the specified file, then loads the specified msp
%           restricted to selMat, which is a multidimensional structure of 
%           frame ids. For all, use datainfo.MspNode(i).Structure
% 
% This function loads the MSOT MSPs from the selected scan folder. 
%
% Return values:
% R          An array containing the MSPed images in 8D
%            1) x
%            2) y
%            3) z (only if 3D detector)
%            4) RUN (repetition of the whole acquisition step, time
%               dimension (timepoints as separate vector)
%            5) z-position (if 2D system with translation stage)
%            6) REPETITION (currently unused)
%            7) Wavelength (see wls vector)
%            8) Individual Images (only if not averaged)
% wls        Wavelength Vector
% zpos       Array of z-stage positions
% ts         Multidimensional array of timestamps in s, same 
%            dimension as R
%%

if ~exist('segm','var')
     % third parameter does not exist, so default it to something
      segm = [];
      toi_exist = 0;
else
    toi_exist = 1;
    mask_toi = segm;
end

datainfo = loadMSOT(msot_file);

% Constant
linewidthTK = 2;
fontsz = 11.5;
map_no = 4;

% this loads the entire experiment
[R, ~, ts, ~] = loadMSOTMsp(msot_file);
msot = squeeze(R);

% obtain the number of spectrally unmixed images 
[d1,d2,d3,d4,d5,n_msp] = size(R);

%% Define ROIs and approve
roi_approve = 0;

% Use ICG signal for ROI identification
prompt = 'Select chromophore for segmentation: ';
str = input(prompt,'s');
cl = str2double(str);
I = squeeze( msot(:,:,:,cl)); 

while roi_approve == 0
    
    % AIF
    ref_anat = I(:,:,50); 
    M=0.5*(max(ref_anat(:)));
    % figure(1); imshow(ref_anat,[0, M], 'Colormap',jet); title('Draw AIF'); h=drawfreehand;
    figure(1); imagesc(ref_anat); colormap(jet); axis image off; title('Draw AIF');
    set(gca, 'ydir', 'normal');
    set(gca, 'xdir', 'reverse');
    h=drawfreehand;
    mask_aif = createMask(h); 

    str = input('Are you satisfied with the AIF ROI? Y/N [Y]: ','s');
    if str == 'Y'
        roi_approve = 1;
    end
end

roi_approve = 0;
while roi_approve == 0

    if toi_exist == 0
        %TOI
        % figure(2); imshow(ref_anat,[0, M], 'Colormap',jet); title('Draw tumor'); h=drawfreehand;
        figure(2); imagesc(ref_anat); colormap(jet); axis image off; title('Draw tumor');
        set(gca, 'ydir', 'normal');
        set(gca, 'xdir', 'reverse');
        h=drawfreehand;
        mask_toi = createMask(h); 

        str = input('Are you satisfied with the Tumor ROI? Y/N [Y]: ','s');
        if str == 'Y'
            roi_approve = 1;
        end
    end

    
end

roi_approve = 0;
while roi_approve == 0

    % Muscle
    % figure(3); imshow(ref_anat,[0, M], 'Colormap', jet); title('Draw muscle'); h=drawfreehand;
    figure(3); imagesc(ref_anat); colormap(jet); axis image off; title('Draw muscle');
    set(gca, 'ydir', 'normal');
    set(gca, 'xdir', 'reverse');
    h=drawfreehand;
    mask_muscle = createMask(h);
    
    str = input('Are you satisfied with the Muscle ROI? Y/N [Y]: ','s');
    if str == 'Y'
        roi_approve = 1;
    end 
   
end

close all;

%Combine masks for saving purposes
if toi_exist == 0
    segm = 3*mask_muscle + 2*mask_aif + mask_toi;
else
    segm = segm + 3*mask_muscle + 2*mask_aif;
end

%% Directory for saving results
cd(saving_DIR);
filename=[datainfo.Name,'_DCE_MSOT_MSP'];
mkdir(filename);
cd(filename);

name = ["deoxy", "oxy", "ICG"];

%% DCE_MSOT fit 
for idx=1:n_msp
    
    I = squeeze( msot(:,:,:,idx) );
    time_min = squeeze( ts(:,:,:,idx) ) / 60;
    [rows,cols, timepoints ] = size(I);

    aif = zeros(1,timepoints);

    for i = 1:timepoints
        S = I(:,:,i);
        aif(i) = mean(S(mask_aif)); 
        % my ROI is too big, I am removing the negative voxels around the AIF;
    end
    
    toi = zeros(1,timepoints);

    for i = 1:timepoints
        S = I(:,:,i);
        toi(i) = mean(S(mask_toi));
    end
    
    muscle = zeros(1, timepoints);
    for i = 1:timepoints
        S = I(:,:,i);
        muscle(i) = mean(S(mask_muscle));
    end
    
    %% Plot average signal for AIF and TOI within ROI
    figure (1);
    hold on;
    plot(time_min, aif, 'bo-', 'LineWidth',linewidthTK);
    plot(time_min, toi, 'ro-', 'LineWidth',linewidthTK);
    plot(time_min, muscle, 'go-', 'LineWidth',linewidthTK);
    hold off;
    xlabel('Time (min)');
    ylabel('Signal ');
    legend ({'AIF','TOI', 'Muscle'},'Location','northeast');
    title(strcat('Signal for',{' '},name(idx)));
    savefig(strcat('Signal for',{' '},name(idx),'.fig'));
    
    
    %% Plot normalized signal and results from fitting
    figure(2);
    hold on;
    plot( time_min, normalize_signal(aif,12)', 'bo-', 'LineWidth',linewidthTK); 
    plot( time_min, normalize_signal(toi,12)', 'ro-', 'LineWidth',linewidthTK);
    plot( time_min, normalize_signal(muscle,12)', 'go-', 'LineWidth',linewidthTK);
    hold off;
    xlabel('Time (min)');
    ylabel('Signal ');
    legend ({'AIF','TOI', 'Muscle'},'Location','northeast');
    title(strcat('Normalized signal for',{' '},name(idx)));
    savefig(strcat('Normalized signal for',{' '},name(idx),'.fig'));
    
    norm_signal_aif = normalize_signal(aif,12)';
    norm_signal_toi = normalize_signal(toi,12)';
    norm_signal_muscle = normalize_signal(muscle,12)';
    
    save('Signal.mat', 'time_min', 'aif', 'toi', 'muscle', 'norm_signal_aif', 'norm_signal_toi', 'norm_signal_muscle');
    
    close all
    
    %% Fitting of average signal in ROI
    [kapp_avg_signal, kep_avg_signal, ~, pars_avg_signal] = fit_dcemsot(toi,aif,time_min);
    [kapp_avg_signal_muscle, kep_avg_signal_muscle,~, pars_avg_signal_muscle] = fit_dcemsot(muscle,aif,time_min);
    
    name_file = strcat('Average_signal_Wavelength_',name(idx),'.mat');
    save(name_file, 'kapp_avg_signal', 'kep_avg_signal','pars_avg_signal'); 
    
    name_file = strcat('Muscle_Average_signal_Wavelength_',name(idx),'.mat');
    save(name_file, 'kapp_avg_signal_muscle', 'kep_avg_signal_muscle', 'pars_avg_signal_muscle'); 
    
    %% Voxel-based fit
    signal = reshape( I, [rows*cols,timepoints]);
    segm = reshape(segm, [rows*cols,1]);
    
    % Allocate space for maps
    kapp = zeros((rows*cols),1);
    kep = zeros((rows*cols),1);
    pars = zeros((rows*cols),map_no);
    
    %% Create Parameter maps using Julio function
    % TOI
    for itr = 1:rows*cols
        if segm (itr)==1
            [kapp(itr,1), kep(itr,1), ~, pars(itr,:)] = fit_dcemsot(double(signal(itr,:)),double(aif), time_min');
        end
    end
    
    kapp_nz = mean(nonzeros(kapp));
    kep_nz = mean(nonzeros(kep));
    kapp = reshape(kapp, [rows,cols]);
    kep = reshape(kep, [rows,cols]);
    pars = reshape(pars, [rows,cols,map_no]);
    segm = reshape(segm,[rows,cols]);
    rsquare_fit = pars(:,:,4);
    kapp_avg = mean(kapp(segm==1));
    kep_avg = mean(kep(segm==1));
    
    name_file = strcat('Maps_Wavelength_',name(idx),'.mat');
    save(name_file, 'kapp_nz', 'kep_nz', 'kapp_avg', 'kep_avg', 'kapp', 'kep', 'rsquare_fit', 'segm','ref_anat'); 
    
    figure;
    ax1 = axes;
    imagesc(ref_anat)
    axis image off; set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
    colormap(ax1, 'gray');
    ax2 = axes;
    imagesc(ax2, kapp, 'alphadata', segm==1);
    axis image off; set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
    colormap(ax2,'jet');
    caxis(ax2, [0 1]);
    ax2.Visible = 'off';
    linkprop([ax1 ax2],'Position');
    title('kapp');
    colorbar;
    savefig('kapp.fig');
    % saveas(gcf, 'kapp.tiff');
    print(gcf, '-dtiff', 'kapp.tiff');
    
    figure;
    ax1 = axes;
    imagesc(ref_anat)
    axis image off; set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
    colormap(ax1, 'gray');
    ax2 = axes;
    imagesc(ax2, kep, 'alphadata', segm==1);
    axis image off; set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
    colormap(ax2,'jet');
    caxis(ax2, [0 1]);
    ax2.Visible = 'off';
    linkprop([ax1 ax2],'Position');
    title('kep');
    colorbar;
    savefig('kep.fig');
    % saveas(gcf, 'kep.tiff'); 
    print(gcf, '-dtiff', 'kep.tiff');
   
    
    % Muscle
    kapp = zeros((rows*cols),1);
    kep = zeros((rows*cols),1);
    pars = zeros((rows*cols),map_no);
    
    for itr = 1:rows*cols
        if segm (itr)==3
            [kapp(itr,1), kep(itr,1), ~, pars(itr,:)] = fit_dcemsot(double(signal(itr,:)),double(aif), time_min');
        end
    end
    
    kapp_nz = mean(nonzeros(kapp));
    kep_nz = mean(nonzeros(kep));
    kapp = reshape(kapp, [rows,cols]);
    kep = reshape(kep, [rows,cols]);
    pars = reshape(pars, [rows,cols,map_no]);
    segm = reshape(segm,[rows,cols]);
    rsquare_fit = pars(:,:,4);
    kapp_avg = mean(kapp(segm==3));
    kep_avg = mean(kep(segm==3));
    
    name_file = strcat('Muscle_Maps_Wavelength_',name(idx),'.mat');
    save(name_file, 'kapp_nz', 'kep_nz', 'kapp_avg', 'kep_avg', 'kapp', 'kep', 'rsquare_fit', 'segm','ref_anat'); 
    
    close all;
end
end

%% Functions
function norm_sig = normalize_signal(mysignal, baseline_points)

norm_sig = mysignal - mean(mysignal(1:baseline_points));
norm_sig = norm_sig / max(mysignal); 
end

function [kapp, kep, X, pars] = fit_dcemsot(TOI,AIF,TIME)

aif_norm = normalize_signal(AIF,12)';
toi_norm = normalize_signal(TOI,12)';

X = zeros(length(aif_norm), 2);
X(:,1) = cumtrapz(TIME, aif_norm);
X(:,2) = -1 .* cumtrapz(TIME, toi_norm);
X(:,3) = ones(size(toi_norm));

y = toi_norm;

options = optimset('Display','Notify');
pars = lsqnonneg(X, y, options);

kapp = pars(1);
kep  = pars(2);

% Calculate predicted curve for the tumor
prediction=X*pars; 

% Calculate the Rsquare for the predicted curve and allocate to pars
pars(4)=rsquare(y,prediction);

end

function [r2, rmse] = rsquare(y,f,varargin)
% Compute coefficient of determination of data fit model and RMSE
%
% [r2 rmse] = rsquare(y,f)
% [r2 rmse] = rsquare(y,f,c)
%
% RSQUARE computes the coefficient of determination (R-square) value from
% actual data Y and model data F. The code uses a general version of 
% R-square, based on comparing the variability of the estimation errors 
% with the variability of the original values. RSQUARE also outputs the
% root mean squared error (RMSE) for the user's convenience.
%
% Note: RSQUARE ignores comparisons involving NaN values.
% 
% INPUTS
%   Y       : Actual data
%   F       : Model fit
%
% OPTION
%   C       : Constant term in model
%             R-square may be a questionable measure of fit when no
%             constant term is included in the model.
%   [DEFAULT] TRUE : Use traditional R-square computation
%            FALSE : Uses alternate R-square computation for model
%                    without constant term [R2 = 1 - NORM(Y-F)/NORM(Y)]
%
% OUTPUT 
%   R2      : Coefficient of determination
%   RMSE    : Root mean squared error
%
% EXAMPLE
%   x = 0:0.1:10;
%   y = 2.*x + 1 + randn(size(x));
%   p = polyfit(x,y,1);
%   f = polyval(p,x);
%   [r2 rmse] = rsquare(y,f);
%   figure; plot(x,y,'b-');
%   hold on; plot(x,f,'r-');
%   title(strcat(['R2 = ' num2str(r2) '; RMSE = ' num2str(rmse)]))
%   
% Jered R Wells
% 11/17/11
% jered [dot] wells [at] duke [dot] edu
%
% v1.2 (02/14/2012)
%
% Thanks to John D'Errico for useful comments and insight which has helped
% to improve this code. His code POLYFITN was consulted in the inclusion of
% the C-option (REF. File ID: #34765).

if isempty(varargin); c = true; 
elseif length(varargin)>1; error 'Too many input arguments';
elseif ~islogical(varargin{1}); error 'C must be logical (TRUE||FALSE)'
else c = varargin{1}; 
end

% Compare inputs
if ~all(size(y)==size(f)); error 'Y and F must be the same size'; end

% Check for NaN
tmp = ~or(isnan(y),isnan(f));
y = y(tmp);
f = f(tmp);

if c; r2 = max(0,1 - sum((y(:)-f(:)).^2)/sum((y(:)-mean(y(:))).^2));
else r2 = 1 - sum((y(:)-f(:)).^2)/sum((y(:)).^2);
    if r2<0
    % http://web.maths.unsw.edu.au/~adelle/Garvan/Assays/GoodnessOfFit.html
        warning('Consider adding a constant term to your model') %#ok<WNTAG>
        r2 = 0;
    end
end

rmse = sqrt(mean((y(:) - f(:)).^2));

end
