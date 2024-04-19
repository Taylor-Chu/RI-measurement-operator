function ROP_param = util_gen_ROP(na, Npb, nTimeSamples, rvtype, ROP_type, Nm)
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

    ROP_param = struct();
    ROP_param.use_ROP = use_ROP;

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

        ROP_param.alpha = alpha;
        ROP_param.beta = beta;
        ROP_param.type = ROP_type;
        ROP_param.rvtype = rvtype;
    end

    if strcmp(ROP_type, 'modul')
        Gamma = randsample([-1,1], Nm*nTimeSamples, true);
        Gamma = reshape(Gamma, [Nm,nTimeSamples]);
        ROP_param.Gamma = Gamma;
    end

end