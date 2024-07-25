function param_ROP = util_gen_ROP(na, Npb, nTimeSamples, rvtype, ROP_type, Nm)
    % Generate random projection vectors
    % na: number of antennas
    % Npb: number of projections per time instant
    % nTimeSamples: number of time instants
    % rvtype: type of random variable ('gaussian' or 'unitary')
    % ROP_type: type of ROP ('none', 'separated' or 'batch')
    % Nm: (optional) number of modulations.

    % flags for using ROPs
    if strcmp(ROP_type, 'separated') || strcmp(ROP_type, 'batch') || strcmp(ROP_type, 'dependent') || strcmp(ROP_type, 'modul')
        use_ROP = true;
    elseif strcmp(ROP_type, 'none')
        use_ROP = false;
    else
        error('ROP_type not recognized')
    end

    param_ROP = struct();
    param_ROP.use_ROP = use_ROP;

    if use_ROP
        if strcmp(rvtype,'gaussian')
            alpha = (randn(na,Npb,nTimeSamples)+1i*randn(na,Npb,nTimeSamples))/sqrt(2);
            beta = (randn(na,Npb,nTimeSamples)+1i*randn(na,Npb,nTimeSamples))/sqrt(2);
        elseif strcmp(rvtype,'unitary')
            alpha = exp(1i*2*pi*rand(na,Npb,nTimeSamples));
            beta = exp(1i*2*pi*rand(na,Npb,nTimeSamples));
        else
            error('Unknown random variable type.');
        end

        param_ROP.alpha = alpha;
        param_ROP.beta = beta;
        param_ROP.type = ROP_type;
        param_ROP.rvtype = rvtype;
    end

    if strcmp(ROP_type, 'modul')
        Gamma = randsample([-1,1], Nm*nTimeSamples, true);
        Gamma = reshape(Gamma, [Nm,nTimeSamples]);
        param_ROP.Gamma = Gamma;
    end

end