function [x] = op_partial_sumst(z, iFt, G, param_ROP)
    % Implement the adjoint of the partial sums operator.
    %
    % Parameters
    % ----------
    % z : array
    %       The modulated ROPs
    % Ft : function handle
    %     Fourier transform operator
    % G : function handle
    %     Visibility interpolation operator
    % param_ROP : struct
    %     Structure containing the parameters for applying ROPs on the measurements
    %
    % Returns
    % -------
    % x  : array
    %      The image
 
    N_zp = size(G,2); % Zero-padded image size
    [Q, Np, B] = size(param_ROP.alpha); % Number of antennas, projections, batches
    Nm = size(param_ROP.Gamma,1); % Number of modulations

    % Adjoint of Bernoulli modulations
    z = reshape(z,[Nm,Np]); % (Nm,Np)

    % useful to subselect Gb
    inds = 1:Q^2*B;
    inds = reshape(inds, [B, Q^2]);

    % open Fx \in N_zp
    Fx = zeros(N_zp,1);

    for b = 1:B

        % compute y_b
        y_b = transpose(z) * param_ROP.Gamma(:,b); % (Np,)

        % Apply R_b*
        alpha_b = param_ROP.alpha(:,:,b);
        beta_b = param_ROP.beta(:,:,b);
        vis_b = ROPt(alpha_b,y_b,beta_b); % (Q^2,)
        vis_b = vis_b/sqrt(B); % (Q^2,)

        % Apply W_b*
        %% This is contained in G_b*

        % Apply G_b*
        inds_b = inds(b,:); % Get Gb 
        G_b = G(inds_b,:); 
        Fx_b = G_b' * vis_b; % (N_zp,)

        % Inject in Fx
        Fx = Fx + Fx_b; % (N_zp,)
    end

    x = iFt(Fx); % Inverse Fourier transform
end