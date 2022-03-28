 function x = pcg_bls_exp_ep(A, C, yi, mi, xi, beta, del, niter, varargin)
% function x = pcg_bls_exp_ep(A, C, yi, mi, xi, beta, del, niter, varargin)
%
% Minimize the cost function with regularizer 4 by using preconditioned conjugate
% gradient with backtracking line search method. The edge-preserving potential
% function is hyperbola function.
% cost(xi) = |yi-Ami.*exp(xi)|_W^2 + beta*|Psi(C exp(xi))|^2
%
% in
%   A   [Nd Np]     system fatrix
%   C   [Nc*Np Np]  finite differencing fatrix
%   yi  [Nd 1]      k-space data
%   mi  [Np 1]      magnitude image
%   xi  [Np 1]      phase image
%   beta    [1]     regularization parameter
%   del [1]     transition point of the potential function
%   niter   [1]     # of subiterations for updating x
%
% option
%   subiter [1]     # of subiterations for line search
%   Pre [1] or [Np Np]  preconditioner
%   thr [1]     stop criterion for line search
%
% out
%   x   [1]     updated phase image
%
% Copyright 2012-06-15, Feng Zhao, The University of Michigan
% Written based on pl_pcg_qs_ls.m in Image Reconstruction Toolbox by
% Jeffrey Fessler, The University of Michigan
% 2013-04-11 modified to enable octave by JF

arg.thr = 0.01;
arg.subiter = 4;
arg.Pre = 1;
arg.btr = 1; % do backtracking or not
arg.mask = true(size(xi));
arg.chat = 0;
arg = vararg_pair(arg, varargin);

oldinprod = 0;

for iter = 1:niter
    if arg.chat
        ticker(mfilename, iter, niter);
    end

    % update xi by roughness regularization
    EX_xi = exp(1i * xi);
%   Am0 = A * diag_sp(mi);
    Am0 = A * Gdiag(mi, 'mask', arg.mask); % jf
    Am0_xn = Am0 * EX_xi;
    CEx = C * EX_xi;
    w1 = wei_func(CEx, del);
    GLn = 2 * real(1i * conj(EX_xi) .* (Am0' * (yi - Am0_xn)));
    GRn = beta * real((-1i) * conj(EX_xi) .* (C' * (CEx .* w1))); % gradient for beta*R(x(n))
    Gn = GLn + GRn; % gradient for Psi(x(n))
    pre_G = arg.Pre * Gn; % preconditioned gradient
    newinprod = Gn' * pre_G;

    % update search direction
    if oldinprod == 0
        dn = -pre_G;
        gamma = 0;
    else
        gamma = newinprod / oldinprod; % Fletcher-Reeves
        dn = -pre_G + gamma * dn;
    end
    oldinprod = newinprod;

    % check if descent direction
    if dn' * Gn > 0
        warning 'wrong direction so resetting';
        printm('<ddir,grad>=%g, |ddir|=%g, |grad|=%g', ...
                dn' * Gn, norm(dn), norm(Gn));
        % reset
        dn = -pre_G;
        oldinprod = 0;
    end

    % optimize step size ak
    ak = 0;
    ak_old = 1e5;
    k = 1;

    % jf added ak==0:
    while ((ak == 0) || (abs(ak - ak_old) / abs(ak) > arg.thr)) && (k <= arg.subiter)
        ak_old = ak;
        EX_xi_dn = exp(1i * (xi + ak * dn));
        Am0_xn_dn = Am0 * EX_xi_dn;
        CEx_dn = C * EX_xi_dn;

        aCEx_dn2 = CEx_dn .* conj(CEx_dn);
        tmp5 = 1 + aCEx_dn2 / del^2;
        tmp6 = sqrt(tmp5);
        w2 = 1 ./ tmp6;
        ww3 = -1 / del^2 .* CEx_dn ./ (tmp6 .* tmp5);

        tmp1 = 1i * conj(EX_xi_dn) .* (Am0' * (yi - Am0_xn_dn));
        tmp2 = (-1i) * conj(EX_xi_dn) .* (C' * (CEx_dn .* w2));
        if ak ~= 0
            GLn = 2 * real(tmp1);
            GRn = beta * real(tmp2); % gradient of beta*R(x(n)+ak*d(n))
        end
        Gf = dn' * (GLn + GRn); % gradient of f(a)

        tmp3 = Am0 * (EX_xi_dn .* dn);
        CLn = 2 * (real(dn' * ((-1i * tmp1) .* dn)) + tmp3' * tmp3); % Hessian of L(x(n)+ak*d(n)) plus some unnecessary imaginery part

        tmp4 = C * (EX_xi_dn .* dn);
        HRn21_dn = 1i * (conj(ww3) .* tmp4);
        HRn22_dn = -1i * (CEx_dn .* conj(tmp4));

        CRn = beta * (real(dn' * (-1i * tmp2 .* dn)) + ...
            real(tmp4' * (w2 .* tmp4)) + ...
            real(HRn21_dn') * real(HRn22_dn)); % curvature of beta*R(x(n)+ak*d(n))

        Cf = CLn + CRn;
        % curvature of f(a)
        Cf = abs(Cf); % take its absolute value
        stp = Gf / Cf;

        if arg.btr == 1 % backtracking
            R1x = sum(del^2 * (tmp6 - 1));
            fa0 = norm(yi - Am0_xn_dn)^2 + beta * R1x;

            stp = 2 * stp;
            while 1
                stp = stp / 2;
                ak1 = ak - stp;
                xi_dn1 = xi + ak1 * dn;
                EX_xi_dn1 = exp(1i * (xi_dn1));
                CEX_dn1 = C * EX_xi_dn1;
                tmp7 = 1 + CEX_dn1 .* conj(CEX_dn1) / del^2;
                R1x = sum(del^2 * (sqrt(tmp7) - 1));

                fa1 = norm(yi - Am0 * EX_xi_dn1)^2 + beta * R1x;
                if (fa1 - fa0) < 1e-8
                    break
                end
                fa0 = fa1;
                printm 'backtracking';
            end
        end

        if Cf == 0
            warning 'found exact solution???  step=0 now!?';
            ak = 0;
        else
            ak = ak - stp;
        end

        k = k + 1;
    end
    if ak < 0
        warning 'downhill step?';
    end

    xi = xi + ak * dn; % update x(n)
end
x = xi;
end

% The weighting function of potential function Psi(|t|)
function w = wei_func(t, del)
w = 1 ./ sqrt(1 + (t .* conj(t)) / del^2);
end
