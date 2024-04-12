function [measop, adjoint_measop] = ops_measop(raw_measop, raw_adjoint_measop, W, Wt)
    % Generate the measurement op and its adjoint from a sampling pattern and
    % user input settings
    % operator (adapted from original code associated with
    % :cite:p:`Fessler2003`).
    %
    % Parameters
    % ----------
    %
    % raw_measop : function handle
    %     Raw measurement operator
    % raw_adjoint_measop : function handle
    %     Adjoint of the raw measurement operator
    % W : function handle
    %     Weighting operator
    % Wt : function handle
    %     Adjoint of the weighting operator
    %
    % Returns
    % -------
    % measop : function handle
    %     Measurement operator
    % adjoint_measop : function handle
    %     Adjoint of the measurement operator
    
    measop = @(x) W(raw_measop(x));
    adjoint_measop = @(y) raw_adjoint_measop(Wt(y));
    
end    