function [tau, noise, gdth_img, vis, param_noise] = util_gen_noise(vis_op, adjoint_vis_op, imSize, nvis, noiselevel, nWimag, param_general, path_uv_data, gdth_img)
    % generate noise realization for the measurements.
    %
    % args:
    %   vis_op: operator computing the visibilities
    %   adjoint_vis_op: adjoint of the visibility operator
    %   nvis: number of visibilities
    %   noiselevel: noise level specification
    %   nWimag: visibility weights
    %   param_general: general parameters
    %   path_uv_data: path to the uv data
    %   gdth_img: ground truth image

    param_noise = struct();
    param_noise.noiselevel = noiselevel;
    param_noise.expo_gdth = true;

    if param_noise.expo_gdth
        % dynamic range of the ground truth image
        log_sigma = rand() * (log10(1e-3) - log10(2e-6)) + log10(2e-6);
        sigma = 10^log_sigma;
        param_noise.targetDynamicRange = 1/sigma;
        if param_general.sigma0 > 0
            % Exponentiation of the ground truth image
            pattern = '(?<=_id_)\d+(?=_dt_)';
            id = regexp(path_uv_data, pattern, 'match');
            seed = str2num(id{1});
            rng(seed, 'twister');
            expo_factor = util_solve_expo_factor(param_general.sigma0, sigma);
            fprintf('\nINFO: target dyanmic range set to %g', param_noise.targetDynamicRange);
            gdth_img = util_expo_im(gdth_img, expo_factor);
        end
    end

    %% Generate the noiseless visibilities
    vis = vis_op(gdth_img);

    %% data noise settings
    switch noiselevel
        case 'drheuristic'
            fprintf("\ngenerate noise (noise level commensurate of the target dynamic range) .. ")

            targetDynamicRange = param_noise.targetDynamicRange;
        
            % include weights in the measurement op.
            measop_1 = @(x) (nWimag.*vis_op(x));
            adjoint_measop_1 = @(x) (adjoint_vis_op(nWimag.*x));
            measopSpectralNorm_1 = op_norm(measop_1, @(y) real(adjoint_measop_1(y)), imSize, 10^-4, 500, 0);

            measop_2 = @(x) ((nWimag.^2) .* vis_op(x));
            adjoint_measop_2 = @(x) (adjoint_vis_op((nWimag.^2).*x));
            measopSpectralNorm_2 = op_norm(measop_2, @(y) real(adjoint_measop_2(y)), imSize, 10^-4, 500, 0);

            % correction factor (=1 if no weights)
            eta_correction = sqrt(measopSpectralNorm_2/measopSpectralNorm_1);

            % noise standard deviation heuristic
            tau  = sqrt(2 * measopSpectralNorm_1) / targetDynamicRange /eta_correction;
            
            % noise realization(mean-0; std-tau)
            noise = tau * (randn(nvis,1) + 1i * randn(nvis,1))./sqrt(2);

            % input signal to noise ratio
            isnr = 20 *log10 (norm(vis)./norm(noise));
            fprintf("\ninfo: random Gaussian noise with input SNR: %.3f db", isnr)

        case 'inputsnr'
            fprintf("\ngenerate noise from input SNR  .. ")

            % user-specified input signal to noise ratio
            param_noise.isnr = 40; % in dB
            isnr = param_noise.isnr; 

            % user-specified input signal to noise ratio
            tau = norm(vis) / (10^(isnr/20)) /sqrt( (nvis + 2*sqrt(nvis)));
            noise = tau * (randn(nvis,1) + 1i * randn(nvis,1))./sqrt(2);
    end

end