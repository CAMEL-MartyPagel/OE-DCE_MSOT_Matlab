function diag_tune_sos(fn,seg,tparams,r_roi,c_array,varargin)

%% first load transducer data
load(tparams);

%% optional parameters
% change channel order
if nargin > 5
    channelorder = varargin{1};
else
    channelorder = [];
end

% invert signal
if nargin > 6
    invert = varargin{2};
else
    invert = 0;
end


%% load data
[sigMat, datainfo] = loadData(fn,'xml',[],0);
[zpositions,timestamps, wavelengths] = get_datainfo(datainfo);

if (~isempty(channelorder))
    sigMat = sigMat(:,channelorder,:);
end
if invert
    sigMat = -sigMat;
end;
%% average dat if necessary
if (size(sigMat,3) > 100)
   numavg = size(sigMat,3) / numel(zpositions); 
   sigMat2 = zeros(size(sigMat,1),size(sigMat,2),numel(zpositions));
   for i = 1:numel(zpositions)
      sigMat2(:,:,i) = mean(sigMat(:,:,(i-1)*numavg+1:i*numavg),3);
   end;
   sigMat = sigMat2;
   clear sigMat2;
end

%%
for ii = 1:numel(datainfo.measurement_desc.projections.msot_frame)
    T(ii) = datainfo.measurement_desc.projections.msot_frame{ii}.actual_temp;
end
T = median(T);
c_form = 1.402385 * 1e3 + 5.038813 * T - 5.799136 * 1e-2 * T^2 + 3.287156 * 1e-4 * T^3 - 1.398845 * 1e-6 * T^4 + 2.787860 * 1e-9 * T^5 ;


%% recon using BP
n = 600;
proj = 300;
limits = [r_sensor-r_roi/1000 r_sensor+r_roi/1000];
filt = [500 7000]*1e3;

BP = zeros(n,n,numel(c_array));
for sl = 1:numel(c_array)
    BP(:,:,sl) = reconWrapper([],sigMat(:,:,seg),[],n,proj,r_sensor,angle_sensor,c_array(sl),'direct',4e7,limits,3,filt,0);
end;


%%
cols = ceil(sqrt(numel(c_array)));
rows = ceil(numel(c_array) / cols);
figure;
for j = 1:numel(c_array)
    subplot(rows,cols,j);
    imagesc(BP(:,:,j));
    axis image;
    colormap gray;
    title([num2str(c_array(j)) ', offset: ' num2str(c_array(j)-c_form) ]);
end