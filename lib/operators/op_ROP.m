function [D,Dt] = op_ROP(param_ROP)
    % Create the operator the applies separated ROPs per batch.
    %
    % Parameters
    % ----------
    % param_ROP : struct
    %     Parameters for the ROP operator.
    %   .alpha : 3D array
    %       The left side projection vectors.
    %   .beta : 3D array
    %       The right side projection vectors.
    %   .type : string
    %       The ROP model. 'separated' or 'batch'.
    %
    % Returns
    % -------
    % D : function handle
    %     Forward operator.
    % Dt : function handle
    %     Adjoint operator.
    
    % extract parameters
    alpha = param_ROP.alpha;
    beta = param_ROP.beta;
    type = param_ROP.type;

    if strcmp(type, 'separated')    
        D = @(x) sep_ROP2(alpha, x, beta);
        Dt = @(y) sep_ROPt(alpha, y, beta);
    elseif strcmp(type, 'batch')
        D = @(x) batch_ROP2(alpha, x, beta);
        Dt = @(y) batch_ROPt(alpha, y, beta);
    elseif strcmp(type, 'dependent')
        D = @(x) dep_ROP(alpha, x, beta);
        Dt = @(y) dep_ROPt(alpha, y, beta);
    elseif strcmp(type, 'modul')
        Gamma = param_ROP.Gamma;
        D = @(x) modul_ROP(alpha, x, beta, Gamma);
        Dt = @(y) modul_ROPt(alpha, y, beta, Gamma);
    else
        error('Unknown ROP type');
    end

end
    