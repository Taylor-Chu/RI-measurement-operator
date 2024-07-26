function param_weighting = util_gen_imaging_weights(param_uv, N, param_weighting)
    %
    % Parameters
    % ----------
    % param_uv : struct
    %   u : double[:]
    %       u coordinate of the data point in the Fourier domain in radians.
    %   v : double[:]
    %       v coordinate of the data points in the Fourier domain in radians.
    % N:  int[2]
    %      image dimension
    % param_weighting : struct
    %     List of parameters to specify weights generation (can be omitted by
    %     default). Fields: `weight_type` with values in {'uniform','briggs', 'none'}, `weight_gridsize`, `weight_robustness` that is the Briggs param.
    %     nW : double[:]
    %       inverse of the noise s.t.d.
    %
    % Returns
    % -------
    % param_weighting : struct
    %   'nWimag' : double[:]
    %      weights inferred from the density of the sampling (uniform/Briggs).
    %
    % author: A. Dabbech, updated [25/05/2024]
    
    %%

    if param_weighting.weighting_on

        if ~isfield(param_weighting, 'weight_type')
            param_weighting.weight_type = 'uniform';
        end
        if ~isfield(param_weighting, 'weight_gridsize')
            param_weighting.weight_gridsize = 2;
        end
        if ~isfield(param_weighting, 'weight_robustness')
            param_weighting.weight_robustness = 0.0;
        end

        imagingBandwidth = param_uv.imagingBandwidth;
        u = param_uv.u.*pi/imagingBandwidth;
        v = param_uv.v.*pi/imagingBandwidth;
        
        % number of visibilities
        nvis = numel(u);
        
        % size of the grid
        N = floor(param_weighting.weight_gridsize*N);
        
        % consider only half of the plane
        u(v < 0) = -u(v < 0);
        v(v < 0) = -v(v < 0);
        
        % grid uv points
        q = floor((u + pi)*N(2)/2/pi);
        p = floor((v + pi)*N(1)/2/pi);
        
        uvInd = sub2ind(N, p, q);
        clear p q;
        
        % Initialize gridded weights matrix with zeros
        gridded_weights = zeros(N);
        
        % inverse of the noise variance
        nW = ones(size(u));
        nW2 = double(nW.^2);
        if isscalar(nW2)
            nW2 = (nW2) .* ones(nvis, 1);
        end
        
        % get gridded weights
        for imeas = 1:nvis
            switch param_weighting.weight_type
                case 'uniform'
                    gridded_weights(uvInd(imeas)) = gridded_weights(uvInd(imeas)) + 1;
                case 'briggs'
                    gridded_weights(uvInd(imeas)) = gridded_weights(uvInd(imeas)) + nW2(imeas);
            end
        
        end
        
        % Apply weighting based on weighting_type
        switch param_weighting.weight_type
            case 'uniform'
                nWimag = 1 ./ sqrt(gridded_weights(uvInd));
            case 'briggs'
                % Compute robust scale factor
                robust_scale = (sum(gridded_weights, "all") / sum(gridded_weights.^2, "all")) * (5 * 10^(-param_weighting.weight_robustness)).^2;
                nWimag = 1 ./ sqrt(1+robust_scale.*gridded_weights(uvInd));
        end
    
    else
        nWimag = ones(length(param_uv.u), 1);
    end

    param_weighting.nWimag = nWimag;
    
end
    