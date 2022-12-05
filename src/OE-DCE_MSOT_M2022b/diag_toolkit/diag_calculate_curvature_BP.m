function [fwhm,fwhm_dist,mip,imax,ymax,norm] = diag_calculate_curvature_BP(fn,numpeaks,mindist,tparams,c_offset,sensors,simfile,outfile,plotres,varargin)
% [fwhm,mip_z] = diag_calculate_curvature_BP(fn,numpeaks,mindist,tparams,c_offset,sensors)
%
% INPUT:
% 1) fn:        Filename of dataset
% 2) numpeaks:  Number of spheres in the phantom
% 3) mindist:   Minimum distance between points
% 4) tparams:   Filename of transducer characteristics (e.g. 256A101)
% 5) c_offset:  Speed of Sound offset using water temperature (or real)
% 6) sensors:   Sensor elements to use for characterisation
% 7) simfile:   Simulation file with specified parameters
% 8) outfile:   Output filename to save results (or [])
% 9) plotres:   Directly plot using diag_plot_curvature
% 10) imasim:   Imasonic simulation file according to surface measurements
 

%% transducer data
load('-mat',tparams);

%% optional parameters
% Imasonic simulation
if nargin > 9
    imasim = varargin{1};
else
    imasim = [];
end

% change channel order
if nargin > 10
    channelorder = varargin{1};
else
    channelorder = [];
end

% invert signal
if nargin > 11
    invert = varargin{2};
else
    invert = 0;
end


%% load data
if (isstruct(fn))
    datainfo = fn;
else
    [datainfo] = loadMSOT(fn);
end
[sigMat] = squeeze(loadMSOTSignals(datainfo));
zpositions = datainfo.ZPositions;
% if (numel(zpositions) > 1 && numel(zpositions) ~= numel(datainfo.ScanFrames))
%     error('Number of Frames and Number of Positions is inconsistent');
% end

if (~isempty(channelorder))
    sigMat = sigMat(:,channelorder,:);
end
if invert
    sigMat = -sigMat;
end;
% %% average dat if necessary
% if (size(sigMat,3) > 200)
%    numavg = size(sigMat,3) / numel(zpositions); 
%    sigMat2 = zeros(size(sigMat,1),size(sigMat,2),numel(zpositions));
%    for i = 1:numel(zpositions)
%       sigMat2(:,:,i) = mean(sigMat(:,:,(i-1)*numavg+1:i*numavg),3);
%    end;
%    sigMat = sigMat2;
%    clear sigMat2;
% end

%% parameters
% recon parameters
n = 300;
limits = [r_sensor-0.02 r_sensor+0.02];
filt = [750 7000]*1e3;

avgwin = 7;         % window for local mip to find normalised map

% for plot labelling
zcenter = floor(size(sigMat,3)./2); % central slice
zres = 0.1;             % resolution along z
ztickdist = 1;          % distance of tickmarks
znum = floor(zcenter*zres/ztickdist);
ztick = zcenter-znum*ztickdist/zres:ztickdist/zres:zcenter+znum*ztickdist/zres;
zticklabel = -znum*ztickdist:ztickdist:znum*ztickdist;


tdist = linspace(limits(1),limits(2),n);
xtickdist = 10;
xticklabel = limits(1)*1e3:xtickdist:limits(2)*1e3;
xtick = round(linspace(1,n,numel(xticklabel)));

%% speed of sound

if c_offset < 100
    T = datainfo.AverageTemperature;
    c = 1.402385 * 1e3 + 5.038813 * T - 5.799136 * 1e-2 * T^2 + 3.287156 * 1e-4 * T^3 - 1.398845 * 1e-6 * T^4 + 2.787860 * 1e-9 * T^5 ;
    % add temperature offset
    c = c+c_offset;   % offset c_form by 22 for phantom
else
    c = c_offset;
end;


%% recon half signals

% create new angles vector that is centered around 270°;
newangles = angle_sensor(sensors);
ccenter = (newangles(1)+newangles(end))./2; % current center angle
newangles = newangles - ccenter +pi;            % should be centered around 0

BP = zeros(n,n,size(sigMat,3));
for sl = 1:size(sigMat,3)
    BP(:,:,sl) = reconWrapper([],sigMat(:,sensors,sl),[],n,numel(sensors)*3,r_sensor,newangles,c,'direct',4e7,limits,3,filt,0);
end;
clear sigMat;


%% MIP in z as well and find peaks
mip = squeeze(max(BP(:,:,5:size(BP,3)-5),[],1));
[mipz z_center] = max(mip,[],2);
% verz rough depth correction
normfreq = linspace(20,1,n)';
mipz = mipz ./ normfreq;

% find extrema and sort by intensity
[ymax imax] = extrema(mipz);
[ymax, I] = sort(ymax,'descend'); imax = imax(I);

% remove double peaks (keep most intense one))
sel = logical(zeros(size(ymax')));
cc = 0;
while nnz(sel) < numpeaks   % loop until all peaks are selected
    cc = cc+1;      if cc > numel(imax) break; end;
    sel(cc) = 1;        % automatically select
    
    % if outside imaging window
    if (imax(cc) > n-avgwin)
        sel(cc) = 0;
    end
    
    for jj = find(sel) % check all previously selected
        % if distance is smaller than mindist
        if ((jj ~= cc) && (abs(diff([imax(jj) imax(cc)])) < mindist))
            sel(cc) = 0;        % then deselect
        end
    end
    
%     if (sel(cc))
%         fprintf('Selecting number %i\n',imax(cc));
%     else
%         fprintf('Rejecting number %i\n',imax(cc));
%     end
end

ymax = ymax(sel); imax = imax(sel);
[imax, I] = sort(imax); ymax = ymax(I);




%% average and condense
clear norm cond;
for i = 1:size(imax)
    cond(i,:) = max(mip(imax(i)-avgwin:imax(i)+avgwin,:),[],1);
    normsph(i,:) = cond(i,:) ./ max(cond(i,:));
end

%% find FWHM from normalised map

loc = normsph >= 0.5;
for i = 1:size(loc,1)
%     fwhm(i) = nnz(loc(i,:));
    fwhm(i) = max(find(loc(i,:)))-min(find(loc(i,:)));
end;
fwhm_dist = tdist(imax);




%% save output
if (~isempty(outfile))
    save('-mat',outfile,...
        'fn','numpeaks','mindist','tparams','c_offset','sensors',...
        'fwhm_dist','fwhm','limits','BP','mip','xtick','xticklabel',...
        'zres','ztick','zticklabel','normsph','ymax','imax','simfile','imasim');
end

%%
if plotres == 1
    diag_plot_curvature(outfile);
end;

return;

%% *** DISPLAY ****

%% plot central image
figure;
subplot(1,2,1);
imagesc(BP_left(:,:,29));
colormap gray;
axis image

subplot(1,2,2);
imagesc(BP_right(:,:,29));
colormap gray;
axis image;

%% look at true MIP in y plane

figure;
imagesc(mip);
% mark peaks
for i = 1:numel(imax)
    line([0 size(mip,1)],[1 1]*imax(i),'Color','white','LineStyle',':');
end


%% plot condensed MIP
figure;
subplot(1,2,1);
imagesc(l_cond);
subplot(1,2,2);
imagesc(r_cond);

%% plot condensed normalised MIP
figure;
subplot(1,2,1);
imagesc(l_norm);
subplot(1,2,2);
imagesc(r_norm);



%% plot profile in x
figure;
plot(lmipz); hold all;
title('profile in x plane');
% plot(rmipz);

for i = 1:numpeaks
    line([l_imax(i) l_imax(i)],[0 max(l_ymax)],'Color','black','LineStyle',':');
end

%% plot FWHM graph

figure;
subplot(1,2,1);
plot(l_tdist(l_imax),l_fwhm./2); hold on;
plot(l_tdist(l_imax),-l_fwhm./2); hold on;

subplot(1,2,2);
plot(r_tdist(r_imax),r_fwhm./2); hold on;
plot(r_tdist(r_imax),-r_fwhm./2); hold on;



