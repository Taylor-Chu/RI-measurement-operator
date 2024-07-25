% Example script to simulate monochromatic RI data from a Fourier sampling pattern
clc; clear ; close all;
fprintf("*** Simulate radio data from a built-in astronomical image ***\n")

%% Setup path
addpath data;
addpath nufft;
addpath lib/operators;
addpath lib/ddes_utils;
addpath lib/utils;

%% ground truth image & settings
% image characteristics
imSize = [512, 512];
% simulations settings
superresolution = 1; % ratio between imaged Fourier bandwidth and sampling bandwidth

% antenna configuration
telescope = 'vlaa';
% total number of snapshots
nTimeSamples = 100;
% obs duration in hours
obsTime = 5;
% obs. frequency in MHz
frequency  = 1e9;
% ROP parameters
use_ROP = true; % use rank-one projections
Npb = 500; % number of projections per time instant
Nm = 100; % number of modulations
ROP_type = 'modul'; % rank-one projected data. ['none', 'dependent', 'separated', 'batch', 'modul']
rvtype = 'unitary'; % or 'gaussian

%% Fourier sampling pattern 
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

%% resolution parameters
resolution_param.superresolution = superresolution; 
% resolution_param.pixelSize = [];

% measurement operator
fprintf("\nbuild NUFFT measurement operator .. ")
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

% %% compute RI normalization factor  (just for info)
dirac = sparse((imSize(1)/2)+1 , (imSize(2)/2)+1 , 1, imSize(1),imSize(2)) ;
psf = real(adjoint_measop(measop(full(dirac))));
ri_normalization = max(psf,[],'all');
figure,imagesc(psf), axis image, colorbar, title('PSF');