function [I_s, W_bda] = util_BDA(path_uv_data)
    %% extract the coordinates
    load(path_uv_data, 'u_ab', 'v_ab'); %, 'w_ab');
    nTimeSamples = size(u_ab, 1);
    na = 64;
    % TODO: take u and v from param_uv instead of reading from file
    u = u_ab(:);
    v = v_ab(:);

    baselines = sqrt(u.^2 + v.^2);
    baseline_min = min(baselines(baselines > 0));

    baseline_edges = [max(baselines)];
    num_pt_avg = [1];
    while true
        baseline_edge_new = baseline_edges(end)/2;
        if baseline_edge_new < baseline_min
            break
        else
            baseline_edges = [baseline_edges baseline_edge_new];
            num_pt_avg = [num_pt_avg num_pt_avg(end) * 2];
        end
    end
    baseline_edges = [baseline_edges 0];
    baseline_edges = flip(baseline_edges);
    num_pt_avg = flip(num_pt_avg);
    num_pt_avg = num_pt_avg(:);
    bins = discretize(baselines, baseline_edges);

    max_num_per_bin = zeros(length(num_pt_avg), 1);
    for i = 1:nTimeSamples
        baselines_b = sqrt(u_ab(i,:,:).^2 + v_ab(i,:,:).^2);
        bins_b = discretize(baselines_b, baseline_edges);
        [C, ~, ic] = unique(bins_b);
        counts = accumarray(ic, 1);
        num_pt_per_bin = ceil(counts ./ num_pt_avg(C));
        max_num_per_bin(C) = max_num_per_bin(C) + num_pt_per_bin;
    end

    M = sum(max_num_per_bin);
    I_s = sparse(M, length(baselines));
    M_count_start = cumsum(max_num_per_bin);
    M_count = zeros(length(max_num_per_bin), 1);
    for i = 1:nTimeSamples
        baselines_b = sqrt(squeeze(u_ab(i,:,:)).^2 + squeeze(v_ab(i,:,:)).^2);
        bins_b = discretize(baselines_b, baseline_edges);
        [C, ~, ic] = unique(bins_b);
        counts = accumarray(ic, 1);
        num_pt_per_bin = ceil(counts ./ num_pt_avg(C));
        counter = 0;
        for j = 1:length(C)
            if num_pt_per_bin(j) == 0
                continue
            else
                bin_b_j_idx = find(bins_b == C(j));
                for k = 1:num_pt_per_bin(j)
                    M_count_start_j = M_count_start(C(j)) + M_count(C(j));
                    k_start = num_pt_avg(C(j)) * (k - 1) + 1;
                    if num_pt_avg(C(j)) * k > length(bin_b_j_idx)
                        k_stop = length(bin_b_j_idx);
                    else
                        k_stop = num_pt_avg(C(j)) * k;
                    end
                    counter = counter + (k_stop - k_start);
                    I_s(M_count_start_j, na^2 * (i - 1) + bin_b_j_idx(k_start:k_stop)) = 1;
                    M_count(C(j)) = M_count(C(j)) + 1;
                end
            end
        end
    end
    W_bda = 1 ./ full(diag(I_s * I_s'));
    W_bda(isinf(W_bda)) = 0;

end