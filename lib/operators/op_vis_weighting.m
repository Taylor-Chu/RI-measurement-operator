function [A, At] = op_vis_weighting(weights)
    % Operator that applies the visibility weighting.
    %
    % Parameters
    % ----------
    % weights : array_like
    %     Array of visibility weights.
    %
    % Returns
    % -------
    %   A : function handle
    %       Function handle that applies the visibility weighting.
    %   At : function handle
    %       Function handle that applies the adjoint visibility weighting.
    
    A = @(x) sqrt(weights) .* x;
    At = @(x) sqrt(weights) .* x;

end    