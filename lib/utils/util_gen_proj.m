function ROP_param = util_gen_proj(na, Npb, nTimeSamples, rvtype)
    % Generate random projection vectors
    % na: number of antennas
    % Npb: number of projections per time instant
    % nTimeSamples: number of time instants
    % rvtype: type of random variable ('gaussian' or 'unitary')

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

end