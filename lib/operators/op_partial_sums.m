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
 
    Q, Np, B = size(param_ROP.alpha); % Number of antennas, projections, batches
    Nm = size(param_ROP.Gamma,1); % Number of modulations

    % Fourier transform of the image
    Fx = Ft(x);

    % useful to subselect Gb
    inds = 1:Q^2*B;
    inds = reshape(inds, [B, Q^2]);

    % open z \in Np x Nm
    z = zeros(Nm,Np);

    for b = 1:B

        % Get Gb 
        inds_b = inds(b,:);
        Gb = G(inds_b,:); 

        % computes the Q^2 visibilities
        vis = Gb* Fx;
        vis = reshape(vis, [Q, Q]);

        if exist('noise', 'var')
            noise_b = noise(inds_b);
            vis = vis + noise_b;
        end

        % Compute y_b to get the ROP
        A = param_ROP.alpha(:,:,b);
        B = param_ROP.beta(:,:,b);
        y_b = ROP(A, vis, B); % (p,1)

        % Compute the Nm modulated versions
        gamma_b = param_ROP.Gamma(:,b); % Nm
        z_b = gamma_b .* transpose(y_b); % (Nm, Np)

        % Inject in z
        z = z + z_b;

    end

end