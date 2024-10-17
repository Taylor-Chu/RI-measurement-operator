function [I_s, W_bda] = util_BDA(path_uv_data)
    %% extract the coordinates
    % path_uv_data = './data/bda_test/uv_2103.mat';
    load(path_uv_data, 'u_ab', 'v_ab'); %, 'w_ab');
    nTimeSamples = size(u_ab, 1);
    na = 64;

    u = u_ab(:);
    v = v_ab(:);

    % compute the baselines of all uv pairs
    baselines = sqrt(u.^2 + v.^2);
    baseline_min = min(baselines(baselines > 0));

    % compute the 'baseline zones' for the BDA, starting from the maximum, decreasing by a factor of 2
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
    num_pt_avg = num_pt_avg(:); % number of points to average over in each baseline zone
    num_pt_per_bin = zeros(length(num_pt_avg), 1); % number of averaged points in each baseline zone

    % reshape u and v into shape (nTimeSamples, na^2) to get the baselines for each antenna pair
    u = reshape(u_ab, [nTimeSamples, na^2]);
    v = reshape(v_ab, [nTimeSamples, na^2]);

    % loop through Q=na^2 to get the number of averaged points in each baseline zone
    for q = 1:na^2
        baseline_q = sqrt(u(:,q).^2 + v(:,q).^2);
        bins_q = discretize(baseline_q, baseline_edges);
        [C, ~, ic] = unique(bins_q);
        counts = accumarray(ic, 1);
        num_pt_per_bin_q = ceil(counts ./ num_pt_avg(C));
        num_pt_per_bin(C) = num_pt_per_bin(C) + num_pt_per_bin_q;
    end

    M = sum(num_pt_per_bin); % total number of averaged points
    % I_s = sparse(M, length(baselines)); % selection matrix of size (M, Q*B^2)
    M_count_start = cumsum(num_pt_per_bin); % starting index of each baseline zone in the selection matrix
    M_count_start = [1; M_count_start(1:end-1)+1];
    M_count = zeros(length(num_pt_per_bin), 1); % counter for the number of averaged points in each baseline zone
    I_s_i = zeros(nTimeSamples * na^2, 1); % row index of the sparse selection matrix
    I_s_j = zeros(nTimeSamples * na^2, 1); % column index of the sparse selection matrix
    W_bda = zeros(M, 1); % initialise the weights for the BDA to avoid computing I_s * I_s'
    I_s_counter = 1;
    max_M = 0;

    % loop through Q=na^2 to fill in the selection matrix
    for q = 1:na^2
        % bin the baselines of the q-th antenna pair
        baseline_q = sqrt(u(:,q).^2 + v(:,q).^2);
        bins_q = discretize(baseline_q, baseline_edges);
        [C, ~, ic] = unique(bins_q); % C = number of baseline zones
        counts = accumarray(ic, 1); % number of baselines in each baseline zone
        num_pt_per_bin_q = ceil(counts ./ num_pt_avg(C)); % number of averaged points in each baseline zone
        % loop through the baseline zones
        for c = 1:length(C)
            if num_pt_per_bin_q(c) == 0
                continue
            else
                bin_q_c_idx = find(bins_q == C(c)); % indices of the baselines in the c-th baseline zone
                % loop through the (to be) averaged points in the c-th baseline zone
                for k = 1:num_pt_per_bin_q(c)
                    M_count_start_k = M_count_start(C(c)) + M_count(C(c)); % starting row index of the k-th averaged point in the c-th baseline zone
                    k_start = num_pt_avg(C(c)) * (k - 1) + 1; % starting column index of the k-th averaged point in the c-th baseline zone
                    if num_pt_avg(C(c)) * k > length(bin_q_c_idx)
                        k_stop = length(bin_q_c_idx);
                    else
                        k_stop = num_pt_avg(C(c)) * k;
                    end
                    k_length = length(bin_q_c_idx(k_start:k_stop)); % number of baselines to be averaged
                    W_bda(M_count_start_k) = 1 / k_length; % assign weight to the points to be averaged
                    I_s_i(I_s_counter:I_s_counter+k_length-1, 1) = M_count_start_k;
                    max_M = max(max_M, M_count_start_k);
                    I_s_j(I_s_counter:I_s_counter+k_length-1, 1) = (q - 1) * nTimeSamples + bin_q_c_idx(k_start:k_stop); % MATLAB matrix is f-contiguous!
                    I_s_counter = I_s_counter + k_length;
                    M_count(C(c)) = M_count(C(c)) + 1;
                end
            end
        end
    end

    I_s = sparse(I_s_i, I_s_j, ones(length(I_s_i), 1), M, nTimeSamples * na^2);
    assert(full(sum(I_s, 'all')) == length(u_ab(:)), "The number of selected visibilities is not equal to the total number of visibilities.");
    clear I_s_i I_s_j baselines u u_ab v v_ab

    % W_bda = 1 ./ full(diag(I_s * I_s'));
    % W_bda(isinf(W_bda)) = 0;

    % for i = 1:nTimeSamples
    %     baselines_b = sqrt(u_ab(i,:,:).^2 + v_ab(i,:,:).^2);
    %     bins_b = discretize(baselines_b, baseline_edges);
    %     [C, ~, ic] = unique(bins_b);
    %     counts = accumarray(ic, 1);
    %     num_pt_per_bin = ceil(counts ./ num_pt_avg(C));
    %     max_num_per_bin(C) = max_num_per_bin(C) + num_pt_per_bin;
    % end

    % M = sum(max_num_per_bin);
    % I_s = sparse(M, length(baselines));
    % M_count_start = cumsum(max_num_per_bin);
    % M_count = zeros(length(max_num_per_bin), 1);
    % for i = 1:nTimeSamples
    %     baselines_b = sqrt(squeeze(u_ab(i,:,:)).^2 + squeeze(v_ab(i,:,:)).^2);
    %     bins_b = discretize(baselines_b, baseline_edges);
    %     [C, ~, ic] = unique(bins_b);
    %     counts = accumarray(ic, 1);
    %     num_pt_per_bin = ceil(counts ./ num_pt_avg(C));
    %     counter = 0;
    %     for j = 1:length(C)
    %         if num_pt_per_bin(j) == 0
    %             continue
    %         else
    %             bin_b_j_idx = find(bins_b == C(j));
    %             for k = 1:num_pt_per_bin(j)
    %                 M_count_start_j = M_count_start(C(j)) + M_count(C(j));
    %                 k_start = num_pt_avg(C(j)) * (k - 1) + 1;
    %                 if num_pt_avg(C(j)) * k > length(bin_b_j_idx)
    %                     k_stop = length(bin_b_j_idx);
    %                 else
    %                     k_stop = num_pt_avg(C(j)) * k;
    %                 end
    %                 counter = counter + (k_stop - k_start);
    %                 I_s(M_count_start_j, na^2 * (i - 1) + bin_b_j_idx(k_start:k_stop)) = 1;
    %                 M_count(C(j)) = M_count(C(j)) + 1;
    %             end
    %         end
    %     end
    % end
    % W_bda = 1 ./ full(diag(I_s * I_s'));
    % W_bda(isinf(W_bda)) = 0;

end