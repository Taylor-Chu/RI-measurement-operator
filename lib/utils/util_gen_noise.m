function [tau, noise] = util_gen_noise(vis_op, adjoint_vis_op, imSize, meas, noise_param, weight_param)
    % generate noise realization for the measurements.
    %
    % args:
    %   vis_op: operator computing the visibilities
    %   adjoint_vis_op: adjoint of the visibility operator
    %   meas: measurements
    %   noise_param: noise level specification
    %   weight_param: parameters associated to visibility weighting

    if nargin < 6
        weighting_on = false;
    else
        weighting_on = weight_param.weighting_on;
        if weighting_on
            nWimag = weight_param.nWimag;
        end
    end

    nmeas = numel(meas);

    switch noise_param.noiselevel
        case 'drheuristic'
            fprintf("\ngenerate noise (noise level commensurate of the target dynamic range) .. ")

            targetDynamicRange = noise_param.targetDynamicRange;
        
            if weighting_on 
                % include weights in the measurement op.
                measop_1 = @(x) (nWimag.*vis_op(x));
                adjoint_measop_1 = @(x) (adjoint_vis_op(nWimag.*x));
                measopSpectralNorm_1 = op_norm(measop_1, @(y) real(adjoint_measop_1(y)), imSize, 10^-4, 500, 0);

                measop_2 = @(x) ((nWimag.^2) .* vis_op(x));
                adjoint_measop_2 = @(x) (adjoint_vis_op((nWimag.^2).*x));
                measopSpectralNorm_2 = op_norm(measop_2, @(y) real(adjoint_measop_2(y)), imSize, 10^-4, 500, 0);

                % correction factor
                eta_correction = sqrt(measopSpectralNorm_2/measopSpectralNorm_1);

                % noise standard deviation heuristic
                tau  = sqrt(2 * measopSpectralNorm_1) / targetDynamicRange /eta_correction;
            else
                % compute measop spectral norm to infer the noise heuristic
                measopSpectralNorm = op_norm(vis_op, @(y) real(adjoint_vis_op(y)), imSize, 10^-4, 500, 0);
                eta_correction = 1;
                % noise standard deviation heuristic
                tau  = sqrt(2 * measopSpectralNorm) / targetDynamicRange ;
            end
            
            % noise realization(mean-0; std-tau)
            noise = tau * (randn(nmeas,1) + 1i * randn(nmeas,1))./sqrt(2);

            % input signal to noise ratio
            isnr = 20 *log10 (norm(meas)./norm(noise));
            fprintf("\ninfo: random Gaussian noise with input SNR: %.3f db", isnr)

        case 'inputsnr'
            fprintf("\ngenerate noise from input SNR  .. ")

            isnr = noise_param.isnr; 

            % user-specified input signal to noise ratio
            tau = norm(meas) / (10^(isnr/20)) /sqrt( (nmeas + 2*sqrt(nmeas)));
            noise = tau * (randn(nmeas,1) + 1i * randn(nmeas,1))./sqrt(2);
    end

end