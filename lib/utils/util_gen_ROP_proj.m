function ROP_proj = util_gen_ROP_proj(na, Npb, nTimeSamples, rvtype)
    % Generate random projection vectors
    % na: number of antennas
    % Npb: number of projections per time instant
    % nTimeSamples: number of time instants
    % rvtype: type of random variable ('gaussian' or 'unitary')
    % ROP_type: type of ROP ('separated' or 'batch')

    if strcmp(rvtype,'gaussian')
        alpha = (randn(na,Npb,nTimeSamples)+1i*randn(na,Npb,nTimeSamples))/sqrt(2);
        beta = (randn(na,Npb,nTimeSamples)+1i*randn(na,Npb,nTimeSamples))/sqrt(2);
    elseif strcmp(rvtype,'unitary')
        alpha = exp(1i*2*pi*rand(na,Npb,nTimeSamples));
        beta = exp(1i*2*pi*rand(na,Npb,nTimeSamples));
    else
        error('Unknown random variable type.');
    end

    ROP_param = struct();
    ROP_param.alpha = alpha;
    ROP_param.beta = beta;
    ROP_param.type = ROP_type;
    ROP_param.rvtype = rvtype;

end