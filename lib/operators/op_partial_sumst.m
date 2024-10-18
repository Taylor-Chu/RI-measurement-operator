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
    Q, Np, B = size(param_ROP.alpha); % Number of antennas, projections, batches
    Nm = size(param_ROP.Gamma,1); % Number of modulations


    % useful to subselect Gb
    inds = 1:Q^2*B;
    inds = reshape(inds, [B, Q^2]);

    % open z \in Np x Nm
    x = zeros(N_zp,1);

    for b = 1:B

        % Get Gb 
        inds_b = inds(b,:);
        Gb = G(inds_b,:); 

        % TODO:
        % IMPLEMENT THE ADJOINT HERE.

        % % computes the Q^2 visibilities
        % vis = Gb* Fx;
        % vis = reshape(vis, [Q, Q]);

        % % Compute y_b to get the ROP
        % A = param_ROP.alpha(:,:,b);
        % B = param_ROP.beta(:,:,b);
        % y_b = ROP(A, vis, B); % (p,1)

        % % Compute the Nm modulated versions
        % gamma_b = param_ROP.Gamma(:,b); % Nm
        % z_b = gamma_b .* transpose(y_b); % (Nm, Np)

        % % Inject in z
        % z =+ z_b;

    end

    x = iFt(Fx); % Inverse Fourier transform

end