function [measop, adjoint_measop] = ops_measop(vis_op, adjoint_vis_op, weight_param, ROP_param)
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
    % weight_param : struct
    %     Structure containing the parameters for visibility weighting
    % ROP_param : struct
    %     Structure containing the parameters for applying ROPs on the measurements
    %
    % Returns
    % -------
    % measop : function handle
    %     Measurement operator
    % adjoint_measop : function handle
    %     Adjoint of the measurement operator
    
    % (optionally) apply visibility weighting
    if weight_param.weighting_on
        w_vis_op = @(x) W(vis_op(x));
        adjoint_w_vis_op = @(y) raw_adjoint_measop(Wt(y));
    else 
        w_vis_op = vis_op;
        adjoint_w_vis_op = adjoint_vis_op;
    end

    % (optionally) apply ROPs
    if ROP_param.use_ROP
        %% compute ROP operator
        [D, Dt] = op_ROP(ROP_param);
    
        measop = @(x) ( D(w_vis_op(x)) ) ; 
        adjoint_measop = @(y) real(adjoint_w_vis_op(Dt(y)));
    else
        measop = w_vis_op;
        adjoint_measop = adjoint_w_vis_op;
    end
    
end    