function [measop, adjoint_measop] = ops_measop(vis_op, adjoint_vis_op, weighting_on, param_ROP)
    % Generate the measurement op and its adjoint from a sampling pattern and
    % user input settings
    % operator (adapted from original code associated with
    % :cite:p:`Fessler2003`).
    %
    % Parameters
    % ----------
    %
    % vis_op : function handle
    %     Visibility operator
    % adjoint_vis_op : function handle
    %     Adjoint of the visibility operator
    % weighting_on : bool
    %     Structure containing the parameters for visibility weighting
    % param_ROP : struct
    %     Structure containing the parameters for applying ROPs on the measurements
    %
    % Returns
    % -------
    % measop : function handle
    %     Measurement operator
    % adjoint_measop : function handle
    %     Adjoint of the measurement operator
    
    % (optionally) apply visibility weighting
    if weighting_on
        w_vis_op = @(x) W(vis_op(x));
        adjoint_w_vis_op = @(vis) adjoint_vis_op(Wt(vis));
    else 
        w_vis_op = vis_op;
        adjoint_w_vis_op = adjoint_vis_op;
    end

    % (optionally) apply ROPs
    if param_ROP.use_ROP
        %% compute ROP operator
        [D, Dt] = op_ROP(param_ROP);
    
        measop = @(x) ( D(w_vis_op(x)) ) ; 
        adjoint_measop = @(y) real(adjoint_w_vis_op(Dt(y)));
    else
        measop = w_vis_op;
        adjoint_measop = adjoint_w_vis_op;
    end
    
end    