function [D,Dt] = op_ROP(ROP_param)
    % Create the operator the applies separated ROPs per batch.
    %
    % Parameters
    % ----------
    % ROP_param : struct
    %     Parameters for the ROP operator.
    %   .Na : int
    %       Number of antennas.
    %   .Npb : int
    %       Number of projections per batch.
    %   .B : int
    %       Number of batches.
    %   .rvtype : string
    %       Type of random variables. Options are 'circ' and 'gaussian'.
    %
    % Returns
    % -------
    % D : function handle
    %     Forward operator.
    % Dt : function handle
    %     Adjoint operator.

    % if ROP_param.rvtype does not exist, set as 'gaussian' by default.
    if ~isfield(ROP_param,'rvtype')
        ROP_param.rvtype = 'gaussian';
    end
    
    % extract parameters
    Na = ROP_param.Na;
    Npb = ROP_param.Npb;
    B = ROP_param.B;

    % generate the random realizations.
    if strcmp(ROP_param.rvtype,'gaussian')
        alpha = (randn(Na,Npb,B)+1i*randn(Na,Npb,B))/sqrt(2);
        beta = (randn(Na,Npb,B)+1i*randn(Na,Npb,B))/sqrt(2);
    elseif strcmp(ROP_param.rvtype,'circ')
        alpha = exp(1i*2*pi*rand(Na,Npb,B));
        beta = exp(1i*2*pi*rand(Na,Npb,B));
    else
        error('Unknown random variable type.');
    end

    % D = @(x) sep_ROP(alpha, x, beta);
    % Dt = @(y) sep_ROPt(alpha, y, beta);
    D = @(x) batch_ROP(alpha, x, beta);
    Dt = @(y) batch_ROPt(alpha, y, beta);

end
    