% Example script to simulate monochromatic RI data from a Fourier sampling pattern
clc; clear ; close all;
fprintf("*** Simulate radio data from a built-in astronomical image ***\n")

%% Setup path
addpath data;
addpath nufft;
addpath lib/operators;
addpath lib/ddes_utils;

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
Npb = 500; % number of projections per time instant
ROP_type = 'separated'; % rank-one projected data. ['none', 'separated', 'batch']
rvtype = 'unitary'; % or 'gaussian

%% flag for using ROPs
if strcmp(ROP_type, 'separated') || strcmp(ROP_type, 'batch') || strcmp(ROP_type, 'dependent')
    use_ROP = true;
elseif strcmp(ROP_type, 'none')
    use_ROP = false;
else
    error('ROP_type not recognized')
end

%% Fourier sampling pattern 
% generate sampling pattern (uv-coverage)
fprintf("\nsimulate Fourier sampling pattern using %s .. ", telescope)
[u, v, w, na] = generate_uv_coverage(frequency, nTimeSamples, obsTime, telescope, use_ROP);

%% generate meas. op & its adjoint
fprintf("\nbuild NUFFT measurement operator .. ")
resolution_param.superresolution = superresolution; 
% resolution_param.pixelSize = [];

% ROP parameters
ROP_param = struct();
if use_ROP
    % generate the random realizations.
    ROP_param = util_gen_ROP(na, Npb, nTimeSamples, rvtype, ROP_type);
end 

% measurement operator
[measop, adjoint_measop] = ops_raw_measop(u,v,w, imSize, resolution_param, ROP_param);

%% perform the adjoint test
measop_vec = @(x) ( measop(reshape(x, imSize)) ); 
adjoint_measop_vec = @(y) reshape(adjoint_measop(y), [prod(imSize), 1]);
measop_shape = struct();
measop_shape.in = [prod(imSize), 1];
if strcmp(ROP_type, 'separated')
    measop_shape.out = [Npb*nTimeSamples,1];
elseif strcmp(ROP_type, 'batch')
    measop_shape.out = [Npb,1];
end
adjoint_test(measop_vec, adjoint_measop_vec, measop_shape);

% %% compute RI normalization factor  (just for info)
dirac = sparse((imSize(1)/2)+1 , (imSize(2)/2)+1 , 1, imSize(1),imSize(2)) ;
psf = real(adjoint_measop(measop(full(dirac))));
ri_normalization = max(psf,[],'all');
figure,imagesc(psf), axis image, colorbar, title('PSF');