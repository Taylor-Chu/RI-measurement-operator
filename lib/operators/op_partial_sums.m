function [z] = op_partial_sums(x, Ft, G, param_ROP, noise)
    % Implement the partial sums view of the forward imaging operator.
    %
    % Parameters
    % ----------
    % x  : array
    %      The input image
    % Ft : function handle
    %     Fourier transform operator
    % G : function handle
    %     Visibility interpolation operator
    % param_ROP : struct
    %     Structure containing the parameters for applying ROPs on the measurements
    % Noise : array
    %     (optional) noise to be added to the visibilities
    %
    % Returns
    % -------
    % z : array
    %       The modulated ROPs
 
    [Q, Np, B] = size(param_ROP.alpha); % Number of antennas, projections, batches
    Nm = size(param_ROP.Gamma,1); % Number of modulations

    % Fourier transform of the image
    Fx = Ft(x);

    % useful to subselect Gb
    inds = 1:Q^2*B;
    inds = reshape(inds, [B, Q^2]);

    % open z \in Nm x Np
    z = zeros(Nm,Np);

    for b = 1:B

        % Get Gb 
        inds_b = inds(b,:);
        G_b = G(inds_b,:); 

        % computes the Q^2 visibilities
        vis_b = G_b * Fx;
        if exist('noise', 'var')
            noise_b = noise(inds_b);
            vis_b = vis_b + noise_b;
        end
        vis_b = reshape(vis_b, [Q, Q]);

        % Compute y_b to get the ROP
        alpha_b = param_ROP.alpha(:,:,b);
        beta_b = param_ROP.beta(:,:,b);
        y_b = ROP(alpha_b, vis_b, beta_b); % (Np,1)

        % Compute the Nm modulated versions
        gamma_b = param_ROP.Gamma(:,b); % (Nm,)
        z_b = gamma_b .* transpose(y_b); % (Nm, Np)

        % Inject in z
        z = z + z_b;

    end

    z = z/sqrt(B); % Normalization

end