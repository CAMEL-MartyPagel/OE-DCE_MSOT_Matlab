function NoiseCompensation( varargin )
% Noise_Compensation  Suppresses parasitic noise in optoacoustic
% signals.
%
%   NOISE_COMPENSATION_ICA()       opens Open File Dialog to choose file
%   NOISE_COMPENSATION_ICA(file)   opens the specified file
% 
% This function loads the signal matrices of a scan and suppresses parasitic
% signals by using independent component analysis.
%
% Required Parameters:
%   <none>
%
% Optional:
%   1: filename    path and name of the msot file to load.
%
% Return Values:
%   <none>
%
% Example:
%   loadMSOT('Scan_1\Scan_1.msot');
%
% Version 1.1 



%% load scan meta data

if numel(varargin) > 0, datainfo = loadMSOT(varargin{1});; else datainfo = loadMSOT(); end;

 
%% perform noise suppression

wbar = waitbar(0,'Performing noise suppression...');

% initialize filters
f_HPF = 50e3;         
fs=40e6;
[b_HPF,a_HPF] = cheby1( 4, .01, 2 * f_HPF/fs * 1.46, 'high' ) ;

if exist([datainfo.FileName,'_save'], 'file') ~= 2
    
    %make a copy of the bin
    copyfile(datainfo.FileName,[datainfo.FileName,'_save'])
    
    %check file size
    fdir_in=dir([datainfo.FileName,'_save']);
    fsize_in=fdir_in.bytes;
    
    %open output and input files
    FIDin = fopen([datainfo.FileName,'_save'],'r');
    FIDout = fopen(datainfo.FileName,'w');
    

    
    % set some parameters
    numproj = datainfo.HWDesc.NumDetectors;
    numsamples = datainfo.MeasurementDesc.RecordLength;
    thresh_corr_coeff=0.85; %threshold for similarity of channels in correlation matrix
    thresh_corr_coeff2=0.7; %threshold for similarity of extracted parasitic noise to the last part of the signal
    truncmin=30;
    sample_compare_min=numsamples-truncmin-200; %first sample for correlation of parasitic noise and signal 
    sample_compare_max=numsamples-truncmin; %last sample for correlation of parasitic noise and signal 
    
    
    %preallocate
    sigMat_out=zeros(numsamples-truncmin,numproj); %init outpur matrix
    ica_corr_val=0; % init correlation value of extracted parasitic noise and the signal

    %loop over frames
    while ~feof(FIDin)
        
        % read scan
        T0=fread(FIDin,[numsamples numproj],'uint16');
        
        %determine offset
        meanvec=mean(T0(100:end,:));
    
        % apply high pass filter
        T = FiltFiltM( b_HPF, a_HPF, T0, 1, 2 );

        % input and output matrices
        sigMat_in=T(truncmin+1:end,:); %truncate input matrix

        % correlation coefficients
        X0=(corrcoef(sigMat_in)); %calculate correlation coefficients to determine similarity between signals in time domain
        X=X0.*(abs(X0)>thresh_corr_coeff).*abs(eye(numproj)-1); % substract trivial self-similarity

        %correction procedure 
        similar_channel_num=zeros(numproj,1); %number of similar channels
        
        %copy input to output signal matrix
        if size(sigMat_in)~=[0,0] %check if input signal matrix is not empty
            sigMat_out=sigMat_in; 
        else
            sigMat_out=zeros(numsamples-truncmin,numproj); %otherwise fill it with zeros
        end
        
        %loop of channels
        if size(sigMat_in)~=[0,0] %check if input signal matrix is not empty
            for channel_num=1:numproj

                if size(X)~=[0,0] %check if correlation matrix is non-zero
                    m=sum(X(channel_num,:)>0); % calculate the number of contributing channels from correlation matrix
                    similar_channel_num(channel_num)=m;
                    clear sigMatTEMP
                    clear sigMatTEMP0
                    sigMatTEMP0=reshape(sigMat_in(repmat(X(channel_num,:),numsamples-truncmin,1)>0),numsamples-truncmin,m);
                    sigMatTEMP=cat(2,sigMat_in(:,channel_num),sigMatTEMP0);

                    if (sum(X(channel_num,:)))>1 %check if similar channels were detected

                        %do ica
                        info = struct;
                        info.factorICA = 1;
                        info.progress = [];
                        [icasig,A,W] = fastica (sigMatTEMP', info, 'lastEig', 1, 'numOfIC', 1,'displayMode', 'off','verbose', 'off','epsilon',0.1);
%                         [icasig,A,W] = fastica (sigMatTEMP', 'lastEig', 1, 'numOfIC', 1,'displayMode', 'off','verbose', 'off','epsilon',0.1);

                        if size(icasig,1)~=0
                            %main component
                            sigcorr=A(1)*icasig';

                            %correlate components with last part of the signal where
                            %no optoacoustic signal is expected
                            ica_corr_mat=corrcoef(cat(2,sigMat_in(sample_compare_min:sample_compare_max,channel_num),icasig(:,sample_compare_min:sample_compare_max)'));

                            %write important correlation matrix elements in a vector
                            ica_corr_val=ica_corr_mat(1,2);

                            %if determined correction exceeds similarity threshold with last part
                            %of the signal, then substract the correction from signal
                            if (abs(ica_corr_val)>thresh_corr_coeff2) 

                                %substract
                                sigMat_out(:,channel_num)=sigMat_in(:,channel_num)-sigcorr;
                            end       
                        end 
                    end
                end

                %set the correlation value to 0
                ica_corr_val=0;

            end
        end
        
        if size(sigMat_in)~=[0,0] %check if input signal matrix is not empty
            fwrite(FIDout,uint16(cat(1,zeros(truncmin,numproj),sigMat_out)+repmat(meanvec,2030,1)),'uint16');
        end
        
        fdir_out=dir([datainfo.FileName]);
        fsize_out=fdir_out.bytes;
        waitbar(fsize_out/fsize_in,wbar,'Performing noise suppression...');
        
    end
    fclose('all');
    close(wbar);
    h = msgbox('Noise suppression completed.','Success');
else
    close(wbar);
    h = msgbox('Operation already performed.', 'Error','error');
end



end
