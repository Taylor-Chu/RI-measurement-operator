function [param_ROP, param_general] = util_gen_ROP(na, nTimeSamples, param_general)
    % Generate random projection vectors
    % na: number of antennas
    % Np: number of projections per time instant
    % nTimeSamples: number of time instants
    % rv_type: type of random variable ('gaussian' or 'unitary')
    % ROP_type: type of ROP ('none', 'separated' or 'batch')
    % Nm: (optional) number of modulations.

    Np = param_general.Np;
    rv_type = param_general.rv_type;
    ROP_type = param_general.ROP_type;
    Nm = param_general.Nm;

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
    param_ROP.Np = Np;
    param_ROP.Nm = Nm;

    if use_ROP
        if strcmp(rv_type,'gaussian')
            alpha = (randn(na,Np,nTimeSamples)+1i*randn(na,Np,nTimeSamples))/sqrt(2);
            beta = (randn(na,Np,nTimeSamples)+1i*randn(na,Np,nTimeSamples))/sqrt(2);
        elseif strcmp(rv_type,'unitary')
            alpha = exp(1i*2*pi*rand(na,Np,nTimeSamples));
            beta = exp(1i*2*pi*rand(na,Np,nTimeSamples));
        else
            error('Unknown random variable type.');
        end

        param_ROP.alpha = alpha;
        param_ROP.beta = beta;
        param_ROP.type = ROP_type;
        param_ROP.rv_type = rv_type;

        % if strcmp(ROP_type, 'separated')
            % param_general.subFolderName = ['separated_ROP', filesep, 'Np_', num2str(Np), param_general.subFolderName]
        % elseif strcmp(ROP_type, 'batch')
            % param_general.subFolderName = ['batch_ROP', filesep, 'Np_', num2str(Np), param_general.subFolderName]
        % elseif strcmp(ROP_type, 'modul')
            % param_general.subFolderName = ['modul_ROP', filesep, 'Np_', num2str(Np), '_Nm_', num2str(Nm), filesep, param_general.subFolderName]
        % end

        if strcmp(ROP_type, 'modul')
            Gamma = randsample([-1,1], Nm*nTimeSamples, true);
            Gamma = reshape(Gamma, [Nm,nTimeSamples]);
            param_ROP.Gamma = Gamma;
        end

    else 
        % param_general.subFolderName = ['no_ROP', filesep, param_general.subFolderName]
    end

end