function img = loadMSOTView(fn,varargin)

if (~exist(fn,'file')),
    error('File does not exist (%s)',fn);
end

% optional parameter: datainfo struct
datainfo = [];
if numel(varargin) > 0,
    datainfo = varargin{1};
end

%% XML parser
try
    patdoc = javaMethod( 'parse', 'com.itheramedical.msotbeans.DataModelImagingSuperSessionDocument$Factory', java.io.File( fn ) ) ;
catch ex
    error( ['Selected XML file could not be loaded. Check well-formedness.\n' ex.message]) ;
end ;


%% View Parameters
ds = patdoc.getDataModelImagingSuperSession.getItems.getDataModelMspImagingSession.getDataSource;
img.SelectedVisibleMainImage = ds.getSelectedVisibleMainImage;

%% get Recon Parameters
rec = patdoc.getDataModelImagingSuperSession.getItems.getDataModelMspImagingSession.getReconPresets.getDataModelReconPreset;
img.recon.n = rec.getSystemReconPresets.getReconstructionPreset.getResolution;
img.OAMPreset.n = rec.getSystemReconPresets.getReconstructionPreset.getResolution;
img.recon.roilow = double(rec.getSystemReconPresets.getReconstructionPreset.getRoiLow);
img.OAMPreset.RoiLow = double(rec.getSystemReconPresets.getReconstructionPreset.getRoiLow);
img.recon.roihigh = double(rec.getSystemReconPresets.getReconstructionPreset.getRoiHigh);
img.OAMPreset.RoiHigh = double(rec.getSystemReconPresets.getReconstructionPreset.getRoiHigh);
img.recon.roicenter = (img.recon.roilow + img.recon.roihigh) / 2;
img.recon.res = (img.recon.roihigh - img.recon.roilow) / img.recon.n;
img.OAMPreset.Projections = rec.getSystemReconPresets.getReconstructionPreset.getProjections;
img.OAMPreset.TimeRes = rec.getSystemReconPresets.getReconstructionPreset.getTimeRes;

img.OAMPreset.FilterLow = rec.getFilterLow;
img.OAMPreset.FilterHigh = rec.getFilterHigh;
img.OAMPreset.DepthCorrection = rec.getDepthCorrection;
img.OAMPreset.ImpulseResponse = rec.getImpulseResponse;
img.OAMPreset.BackgroundAbsorption = double(rec.getBackgroundAbsorption);
img.OAMPreset.BackgroundOxygenation = double(rec.getBackgroundOxygenation);
img.OAMPreset.UserSoundTrim = double(rec.getUserSoundTrim);

% MSP
mp = patdoc.getDataModelImagingSuperSession.getItems.getDataModelMspImagingSession.getMspPresets.getDataModelMspPreset;
img.OAMPreset.MSP.DiscardNegatives = mp.getDiscardNegativeValues;
img.OAMPreset.MSP.Method = char(mp.getMethod);
img.OAMPreset.MSP.BGWavelength = mp.getBgWavelength.getIntValue;
img.OAMPreset.MSP.Spectra = cell(mp.getUserSelectedSpectra.getStringArray);

% Settings
sett = patdoc.getDataModelImagingSuperSession.getItems.getDataModelMspImagingSession.getSettings;
img.OAMPreset.ViewSettings.UltrasoundMinimumScaling = double(sett.getUltrasoundMinimumScaling);
img.OAMPreset.ViewSettings.UltrasoundMaximumScaling = double(sett.getUltrasoundMaximumScaling);
img.OAMPreset.ViewSettings.BackgroundMinimumScaling = double(sett.getBackgroundMinimumScaling);
img.OAMPreset.ViewSettings.BackgroundMaximumScaling = double(sett.getBackgroundMaximumScaling);
img.OAMPreset.ViewSettings.ForegroundMinimumScaling = double(sett.getForegroundMinimumScaling);
img.OAMPreset.ViewSettings.ForegroundMaximumScaling = double(sett.getForegroundMaximumScaling);
img.OAMPreset.ViewSettings.Visible3DGridPlanesTypes = char(sett.getVisible3DGridPlanesTypes);

%% get Image Tags
if ~isempty(datainfo),
    framets = [datainfo.ScanFrames.Timestamp];
    framets = framets - framets(1);
end

tagarr = patdoc.getDataModelImagingSuperSession.getItems.getDataModelMspImagingSession.getDataSource.getImageTaggingList.getDataModelImageTaggingArray;
img.Tag = [];
for t = 1:numel(tagarr),
   img.Tag(t).Timestamp = double(tagarr(t).getTimestamp)*1e-7;
   if ~isempty(datainfo),
       [~,frameind] = min(abs(framets - img.Tag(t).Timestamp));
       img.Tag(t).Repetition = datainfo.ScanFrames(frameind).Repetition;
   else
       img.Tag(t).Repetition = nan;
   end
end

%% get ROI Regions
sett = patdoc.getDataModelImagingSuperSession.getItems.getDataModelMspImagingSession.getSettings;
roiarr = sett.getDrawingRois2D.getRegions.getDataModel3DRegionArray;

% figure,ih = image([img.recon.roilow img.recon.roihigh],[img.recon.roilow img.recon.roihigh],zeros(332,332));axis image;hold on;

roi = cell(numel(roiarr),1);
roidrop = false(numel(roiarr),1);
for r = 1:numel(roiarr),
    roi{r}.name = { char(roiarr(r).getName) };
    col = roiarr(r).getColor.getColorProperty;
    roi{r}.color = [ double(col.getScR)  double(col.getScG)  double(col.getScB) ];
    base = roiarr(r).getRegion3D.getDataModel3DRegionBase;
    roi{r}.class = lower(char(base.getShape));
    roipoints = base.getRegionPoints.getPointArray;
    
    pos = zeros(numel(roipoints),2);
    for ii = 1:numel(roipoints)
%         pos(ii,1) = img.recon.roihigh - img.recon.roicenter - double(roipoints(ii).getX)*img.recon.res;
        pos(ii,1) = img.recon.roilow - img.recon.roicenter + double(roipoints(ii).getX)*img.recon.res;
        pos(ii,2) = img.recon.roilow - img.recon.roicenter + double(roipoints(ii).getY)*img.recon.res;
%         fprintf('%s (%s): (%.1f,%.1f)\n',roi{r}.name,roi{r}.class,pos(ii,1),pos(ii,2));
    end
%     plot(pos(:,1),pos(:,2),'*','Color',roi{r}.color); hold all;
    
    % convert positions
    switch (roi{r}.class)
        case 'rectangle'
            roi{r}.class = 'imrect';
            roi{r}.pos(1) = min(pos(:,1));
            roi{r}.pos(2) = min(pos(:,2));
            roi{r}.pos(3) = abs(diff(pos(:,1)));
            roi{r}.pos(4) = abs(diff(pos(:,2)));
%             ir = imrect(gca,roi{r}.pos);
%             setColor(ir,roi{r}.color);
        case 'circle'
            roi{r}.class = 'imellipse';
            radius = sqrt((pos(2,1)-pos(1,1))^2 + (pos(2,2)-pos(1,2))^2);
            roi{r}.pos(1) = pos(1,1)-radius;
            roi{r}.pos(2) = pos(1,2)-radius;
            roi{r}.pos(3) = 2*radius;
            roi{r}.pos(4) = 2*radius;
%             ie = imellipse(gca,roi{r}.pos);
%             setColor(ie,roi{r}.color);
        case 'ellipse'
            warning('Importing ellipses is currently not supported');
            roidrop(r) = true;
%             
%              a = sqrt((pos(1,1)-pos(2,1))^2+(pos(1,2)-pos(2,2))^2);
%              b = sqrt((pos(3,1)-pos(2,1))^2+(pos(3,2)-pos(2,2))^2);
%              t = linspace(0,2*pi);
%              X = a*cos(t);
%              Y = b*sin(t);
%              w = atan2(pos(1,2)-pos(2,2),pos(1,1)-pos(2,1));
%              x = pos(2,1) + X*cos(w) - Y*sin(w);
%              y = pos(2,2) + X*sin(w) + Y*cos(w);
%              roi{r}.pos(1) = min(x);
%              roi{r}.pos(2) = min(y);
%              roi{r}.pos(3) = max(x)-min(x);
%              roi{r}.pos(4) = max(y)-min(y);
% %              imellipse(gca,roi{r}.pos);
%             plot(x,y,'.');
%              roi{r}.pos
        case 'polygon'
            roi{r}.class = 'impoly';
            roi{r}.pos = pos;
%             ip = impoly(gca,pos);
%             setColor(ip,roi{r}.color);
        otherwise
            warning('Unsupported ROI Class: %s - ignored',roi{r}.class);
    end
    
end


img.roi = roi(~roidrop);