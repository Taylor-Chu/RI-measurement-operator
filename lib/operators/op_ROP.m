function [D,Dt] = op_ROP(ROP_proj)
    % Create the operator the applies separated ROPs per batch.
    %
    % Parameters
    % ----------
    % ROP_proj : struct
    %     Parameters for the ROP operator.
    %   .alpha : 3D array
    %       The left side projection vectors.
    %   .beta : 3D array
    %       The right side projection vectors.
    %
    % Returns
    % -------
    % D : function handle
    %     Forward operator.
    % Dt : function handle
    %     Adjoint operator.
    
    % extract parameters
    alpha = ROP_param.alpha;
    beta = ROP_param.beta;
    type = ROP_param.type;

    if strcmp(type, 'separated')    
        D = @(x) sep_ROP(alpha, x, beta);
        Dt = @(y) sep_ROPt(alpha, y, beta);
    elseif strcmp(type, 'batch')
        D = @(x) batch_ROP(alpha, x, beta);
        Dt = @(y) batch_ROPt(alpha, y, beta);
    end

end
    