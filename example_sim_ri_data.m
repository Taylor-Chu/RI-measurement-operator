% Example script to simulate RI data, using a toy Fourier sampling pattern

% clc; clear ; close all;
fprintf("*** Simulate toy radio data from a built-in astronomical image ***\n")

%% Setup paths
addpath data;
addpath nufft;
addpath lib/operators;
addpath lib/utils;
addpath lib/ddes_utils;

%% simulation setting: realistic / toy
simtype = 'realistic'; % possible values: `realistic` ; `toy`
noiselevel = 'drheuristic'; % possible values: `drheuristic` ; `inputsnr`
superresolution = 1.5; % ratio between imaged Fourier bandwidth and sampling bandwidth

switch simtype
    case 'realistic'
        myuvwdatafile = 'tests/test.mat';
        frequency = load(myuvwdatafile,'frequency').frequency;
        % data weighting enabled (for imaging e.g. Briggs)  
        weighting_on = true; 

    case 'toy'
        % antenna configuration
        telescope = 'vlaa';
        % total number of snapshots
        nTimeSamples = 100;
        % obs duration in hours
        obsTime = 4;
        % obs. frequency in MHz
        frequency  = 1e9;
        % data weighting enabled (for imaging e.g. Briggs) 
        weighting_on = false; 
end

%% ground truth image 
fprintf("\nread ground truth image  .. ")
% built-in matlab image
gdthim = imread('ngc6543a.jpg') ; 
% crop region of interest
gdthim = double(gdthim(49:560, 27:538)); 
% normalize image (peak = 1)
gdthim = gdthim./max(gdthim,[],'all'); 
% characteristics
imSize = size(gdthim);
% display
figure(); imagesc(gdthim); colorbar; title('ground truth image'); axis image;  axis off;

%% data noise settings
noise_param = struct();
noise_param.noiselevel = noiselevel;
switch noiselevel
    case 'drheuristic'
        % dynamic range of the ground truth image
        noise_param.targetDynamicRange = 255; 
    case 'inputsnr'
         % user-specified input signal to noise ratio
        noise_param.isnr = 40; % in dB
end

%% Parameters for visibility weighting
weight_param = struct();
weight_param.weighting_on = weighting_on;

%% Fourier sampling pattern
switch simtype
    case 'realistic'
        fprintf("\nload Fourier sampling pattern .. ")
        uvwdata = load(myuvwdatafile,'u','v','w');
        umeter =  uvwdata.u;
        vmeter =  uvwdata.v;
        wmeter =  uvwdata.w;
        clear uvwdata
        speedOfLight = 299792458;
        u = umeter ./ (speedOfLight/frequency) ;
        v = vmeter ./ (speedOfLight/frequency) ;
        w = wmeter ./ (speedOfLight/frequency) ;
        
        % for info 
        try nominalPixelSize = double(load(myuvwdatafile,'nominal_pixelsize').nominal_pixelsize);
        end

        % imaging weights if available
        if weighting_on
            try nWimag = double(load(myuvwdatafile,'nWimag').nWimag);
                weight_param.nWimag = nWimag(:);
                fprintf("\ninfo: imaging weights are found")
            catch
                weighting_on = false;
                weight_param.weighting_on = weighting_on;
            end
        end

    case 'toy'       
        % generate sampling pattern (uv-coverage)
        fprintf("\nsimulate Fourier sampling pattern using %s .. ", telescope)
        [u, v, w, na] = generate_uv_coverage(frequency, nTimeSamples, obsTime, telescope);
end

% maximum projected baseline (just for info)
maxProjBaseline  = sqrt(max(u.^2+v.^2));

%% generate meas. op & its adjoint
fprintf("\nbuild NUFFT measurement operator .. ")
resolution_param.superresolution = superresolution; 
% resolution_param.pixelSize = nominalPixelSize/superresolution; 
[raw_measop, adjoint_raw_measop] = ops_raw_measop(u,v, w, imSize, resolution_param);
 
%% model clean visibilities 
fprintf("\nsimulate model visibilities .. ")
vis = raw_measop(gdthim);

%number of data points
nmeas = numel(vis);

%% model data

% noise vector
[tau, noise] = util_gen_noise(raw_measop, adjoint_raw_measop, imSize, vis, noise_param, weight_param);

% data
fprintf("\nsimulate data  .. ")
y = vis + noise;

%% Eventually switch visibility weighting on
if weighting_on
    [W, Wt] = op_vis_weighting(weight_param.nWimag);
    [measop, adjoint_measop] = ops_measop(raw_measop, adjoint_raw_measop, W, Wt);
end

%% back-projected data
fprintf("\nget (non-normalised) back-projected data  .. ")
title_ = 'dirty image';
if weighting_on
    dirty = real( adjoint_measop(nWimag.*y) );
    title_ = [title_ , ' (weights applied)'];
else
    dirty = real( adjoint_raw_measop(y) );
end
figure(); imagesc(dirty); colorbar; title(title_); axis image; axis off;

%% generate input data file for uSARA/AIRI/R2D2 imager  (just for info)
% whitening vector
fprintf("\nsave data file  .. ")
mkdir 'results'
matfilename = "results/ngc6543a_data.mat" ;
dirtyfilename = "results/ngc6543a_dirty.fits" ; 

% save mat file
nW = tau *ones(nmeas,1);
save(matfilename, "y", "nW", "u", "v","w","maxProjBaseline","frequency",'-v7.3')
% add imaging weights
if weighting_on
    save(matfilename,"nWimag",'-append')
end

% save (non-normalised) dirty image
fitswrite(dirty, dirtyfilename)


fprintf('\nDone.')
