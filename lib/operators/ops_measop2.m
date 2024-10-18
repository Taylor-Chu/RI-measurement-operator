function [measop, adjoint_measop, y] = ops_measop2(x, G, Ft, IFt, param_weighting, tau, param_ROP)
    % Generate the measurement op and its adjoint from a sampling pattern and
    % user input settings
    % operator (adapted from original code associated with
    % :cite:p:`Fessler2003`).
    %
    % Parameters
    % ----------
    %
    % G : function handle
    %     Visibility interpolation operator
    % Ft : function handle
    %     Fourier transform operator
    % IFt : function handle
    %     Inverse Fourier transform operator
    % param_weighting : struct
    %     weighting_on : bool
    %         Structure containing the parameters for visibility weighting
    %     nWimag : double[:]
    %         weights inferred from the density of the sampling (uniform/Briggs)
    % tau : float
    %     Noise level
    % param_ROP : struct
    %     Structure containing the parameters for applying ROPs on the measurements
    %
    % Returns
    % -------
    % measop : function handle
    %     Measurement operator
    % adjoint_measop : function handle
    %     Adjoint of the measurement operator
    % y : double[:]
    %     Measurement vector
    
    nWimag = param_weighting.nWimag;
    weighting_on = param_weighting.weighting_on;
    precompute = param_ROP.precompute;

    % (optionally) apply visibility weighting
    if weighting_on
        % nW = (1 / tau) * ones(na^2*nTimeSamples,1);
        nW = (1 / tau) * nWimag;
        vis = nW.*vis;

        % Integrate the weighting into the interpolation matrix G.
        G = nW .* G;
    end

    measop = @(x) op_partial_sums(x, Ft, G, param_ROP);

    adjoint_measop = @(z) op_partial_sumst(z, iFt, G, param_ROP);

    y = op_partial_sums(x, Ft, G, param_ROP, noise); % Compute the modulated ROPs from the groundtruth image

    % measop = @(x) ( G * Ft(x) ) ; 
    % adjoint_measop = @(y) real(IFt(G' * y));
    
    % % (optionally) apply ROPs
    % if param_ROP.use_ROP

    %     [D, Dt] = op_ROP(param_ROP);

    %     measop = @(x) ( D(measop(x)) );
    %     adjoint_measop = @(y) adjoint_measop(Dt(y));
    %     y = D(vis);

    %     % Precompute forward operator
    %     if strcmp(param_ROP.type, 'modul')

    %         if precompute
    %             [DG, adjoint_DG] = util_precompute_op(G, param_ROP);

    %             measop = @(x) ( DG(Ft(x)) );
    %             adjoint_measop = @(y) real(IFt(adjoint_DG(y)));
    %         end

    %     end
    % else
    %     y = vis;
    % end
    
end    