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
if strcmp(ROP_type, 'separated') || strcmp(ROP_type, 'batch')
    use_ROP = true;
elseif strcmp(ROP_type, 'none')
    use_ROP = false;
else
    error('ROP_type not recognized')
end

%% Fourier sampling pattern 
% generate sampling pattern (uv-coverage)
fprintf("\nsimulate Fourier sampling pattern using %s .. ", telescope)
[umeter, vmeter, wmeter, na] = generate_uv_coverage(nTimeSamples, obsTime, telescope);

% convert in units of the wavelength
speedOfLight = 299792458;
u = umeter ./ (speedOfLight/frequency) ;
v = vmeter ./ (speedOfLight/frequency) ;
w = wmeter ./ (speedOfLight/frequency) ;
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

% %% compute RI normalization factor  (just for info)
dirac = sparse((imSize(1)/2)+1 , (imSize(2)/2)+1 , 1, imSize(1),imSize(2)) ;
psf = real(adjoint_measop(measop(full(dirac))));
ri_normalization = max(psf,[],'all');
figure,imagesc(psf), axis image, colorbar, title('PSF');