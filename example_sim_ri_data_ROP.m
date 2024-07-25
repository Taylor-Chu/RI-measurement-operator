% Example script to simulate RI data, using a toy Fourier sampling pattern

% clc; clear ; close all;
fprintf("*** Simulate toy rank-one projected radio data from a built-in astronomical image ***\n")

%% Setup paths
addpath data;
addpath nufft;
addpath lib/operators;
addpath lib/operators/ROP/;
addpath lib/utils;
addpath lib/ddes_utils;

%% simulation setting: realistic / toy
noiselevel = 'drheuristic'; % possible values: `drheuristic` ; `inputsnr`
superresolution = 1.5; % ratio between imaged Fourier bandwidth and sampling bandwidth

% antenna configuration
telescope = 'vlaa';
% total number of snapshots
nTimeSamples = 100; 
% obs duration in hours
obsTime = 4;
% obs. frequency in MHz
frequency  = 1e9;
% ROP parameters
use_ROP = true; % use rank-one projections
Npb = 500; % number of projections per time instant
Nm = 100; % number of modulations
ROP_type = 'modul'; % rank-one projected data. ['none', 'dependent', 'separated', 'batch', 'modul']
rvtype = 'unitary'; % or 'gaussian

%% ground truth image 
fprintf("\nread ground truth image  .. ")
% built-in matlab image
gdthim = imread('ngc6543a.jpg') ; 
% crop region of interest
gdthim = double(gdthim(49:560, 27:538)); 
% normalize image (peak = 1)
gdthim = gdthim./max(gdthim,[],'all'); 

imSize = size(gdthim); % characteristics

figure(1), imagesc(gdthim), colorbar, title ('ground truth image'), axis image,  axis off,

% generate sampling pattern (uv-coverage)
fprintf("\nsimulate Fourier sampling pattern using %s .. ", telescope)
[u, v, w, na] = generate_uv_coverage(frequency, nTimeSamples, obsTime, telescope, use_ROP);

% figure(); plot(u, v, 'o'); title('uv-coverage'); axis equal; grid on;

%% uv-coverage data
param_uv = struct();
param_uv.u = u;
param_uv.v = v;
param_uv.w = w;
param_uv.na = na;
param_uv.nTimeSamples = nTimeSamples;

%% Set ROP parameters
param_ROP = util_gen_ROP(na, Npb, nTimeSamples, rvtype, ROP_type, Nm);

%% data noise settings
param_noise = struct();
param_noise.noiselevel = noiselevel;
switch noiselevel
    case 'drheuristic'
        % dynamic range of the ground truth image
        param_noise.targetDynamicRange = 255; 
    case 'inputsnr'
         % user-specified input signal to noise ratio
        param_noise.isnr = 40; % in dB
end

% maximum projected baseline (just for info)
maxProjBaseline  = sqrt(max(u.^2+v.^2));

%% generate meas. op & its adjoint
fprintf("\nbuild NUFFT measurement operator .. ")
resolution_param.superresolution = superresolution; 
% resolution_param.pixelSize = nominalPixelSize/superresolution; 

% measurement operator
[measop, adjoint_measop] = ops_raw_measop(param_uv, imSize, resolution_param, param_ROP);
 
% %% perform the adjoint test
% measop_vec = @(x) ( measop(reshape(x, imSize)) ); 
% adjoint_measop_vec = @(y) reshape(adjoint_measop(y), [prod(imSize), 1]);
% measop_shape = struct();
% measop_shape.in = [prod(imSize), 1];
% if strcmp(ROP_type, 'separated')
%     measop_shape.out = [Npb*nTimeSamples,1];
% elseif strcmp(ROP_type, 'batch')
%     measop_shape.out = [Npb,1];
% elseif strcmp(ROP_type, 'dependent')
%     measop_shape.out = [Npb^2*nTimeSamples,1];
% elseif strcmp(ROP_type, 'modul')
%     measop_shape.out = [Npb*Nm,1];
% end
% adjoint_test(measop_vec, adjoint_measop_vec, measop_shape);

%% model clean measurements
fprintf("\nsimulate model measurements .. ")
meas = measop(gdthim);

nmeas = numel(meas); %number of data points

%% model data

% noise vector
[tau, noise] = util_gen_noise(measop, adjoint_measop, imSize, meas, param_noise);

% data
fprintf("\nsimulate noisy data  .. ")
y = meas + noise;

%% back-projected data
fprintf("\nget (non-normalised) back-projected data  .. ")
dirty = real( adjoint_measop(y) );
dirty = reshape(dirty, imSize);
figure(2), imagesc(dirty), colorbar, title ('dirty image'), axis image,   axis off,

%% generate input data file for uSARA/AIRI/R2D2 imager  (just for info)
% whitening vector
fprintf("\nsave data file  .. ")
mkdir 'results'
if use_ROP
    ROP_text = ['_ROP_', ROP_type, '_', rvtype, num2str(Npb)];
    matfilename = ['results/ngc6543a_data', ROP_text, '.mat'];
    dirtyfilename = ['results/ngc6543a_dirty', ROP_text,'.fits'] ; 
else 
    matfilename = ['results/ngc6543a_data.mat'];
    dirtyfilename = ['results/ngc6543a_dirty.fits'] ; 
end
gtfilename = "results/ngc6543a_gt.fits" ;

% save mat file
nW = tau * ones(na^2*nTimeSamples,1);
save(matfilename, "y", "nW", "u", "v","w","maxProjBaseline","frequency", "param_ROP", '-v7.3')

% save (non-normalised) dirty image
fitswrite(dirty, dirtyfilename)

% save ground truth image
fitswrite(gdthim, gtfilename)


fprintf('\nDone.')
