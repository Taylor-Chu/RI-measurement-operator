[I_s, W_bda] = util_BDA(path_uv_data)
    %% extract the coordinates
    load(path_uv_data, 'u_ab', 'v_ab'); %, 'w_ab');
    u = u_ab(:);
    v = v_ab(:);
    % w = w_ab(:);

    % convert unit to km
    u = u / 1000;
    v = v / 1000;
    % w = w / 1000;

    baselines = sqrt(u.^2 + v.^2);
    da_baselines = [80, 40, 30, 20, 15, 10, 7.5, 5, 3.75, 2.5, 1.875, 1.25, 0.9375, 0.625, 0.5625, 0.375, 0.28125, 0];
    bda_baselines = flip(bda_baselines);
    bda_avg = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512];
    bda_avg = flip(bda_avg);

    bins = discretize(baselines, bda_baselines);

    num_bins = length(unique(bins));

    I_s = zeros([num_bins, length(baselines)]);
    for i = 1:length(bins)
        I_s(bins(i), i) = 1;
    end

    W_bda = 1 ./ (I_s * I_s');
    W_bda(isinf(W_bda)) = 0;

end