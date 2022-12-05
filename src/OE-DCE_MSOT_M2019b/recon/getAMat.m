function [A_mat Amatfile] = getAMat(par,varargin)
% GETAMAT  Loads or creates model matrix for model based reconstruction
%   [A_mat Amat_file] = getAmat(par [, angle_sensor, t]);
%
%   Parameters:
%   - par: Parameter struct containing all important properties for naming
%     * Amat_dir: Directory containing retained model matrices
%                 (will not save if directory is undefined or not existing)
%     * n:        Number of pixels
%     * r:        Radius of transducer
%     * c:        Speed of sound
%     * proj:     Number of projections (interpolated if > nb of elements)
%     * angles:   Number of interpolation angles (model based)
%     * timeres:  Time resolution (model based)
%     * depth:    Decomposition depth (wavelet)
%     * Athres:   Threshold for model matrix elements (wavelet)
%     * Sthres:   Threshold of singular values (wavelet)
%   - angle_sensor: sensor coordinates (optional, req for matrix calc)
%   - t:            time vector (optional, req for matrix calc)
%
%   Return values:
%   - A_mat:      Model matrix (projections*timepoints x pixel^2)
%   - Amat_file:  Filename of model matrix (empty if not saved)
    
    A_mat = [];
    Amatfile = [];
    if (isfield(par,'Amat_dir') && exist(par.Amat_dir,'dir'))
        Amatfile = [par.Amat_dir ...
            '\Amat_roi' num2str(par.roi*1e3,'%02d') ...
            '_n' num2str(par.n) ...
            '_r' num2str(par.r_sensor*1e3,'%.2f') ...
            '_c' num2str(round(par.c)) ...
            '_proj' num2str(par.proj) ...
            '_angles' num2str(par.n_angles) ...
            '_timeres' num2str(par.timeres) ...
            '.mat'];

        % if File already exists, load and return
        if (exist(Amatfile,'file') && isempty(par.limitSensors))
            fprintf('Loading A_mat (%s)...\n',Amatfile);
            load(Amatfile);
            return;
        end
    end
        
    % if file does not yet exist (or should never exist)
    if numel(varargin) < 2
        error('Please supply angle_sensor and t for matric calculation');
        return;
    end
    angle_sensor = varargin{1};
    t = varargin{2};

    fprintf('Calculating and saving A_mat (%s)...\n',Amatfile);
    A_mat = computeDiscretizedMatrix(par.c,par.n,par.roi,t,par.r_sensor,angle_sensor,par.n_angles)';
        

    % if filename exists, save file (otherwise recreate each time)
    if (~isempty(Amatfile))
        save(Amatfile,'-v7.3','A_mat');
    else
        warning('Not saving model matrix, please set Amat_dir');
    end
    