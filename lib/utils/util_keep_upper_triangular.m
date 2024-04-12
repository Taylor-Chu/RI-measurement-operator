function u2 = util_keep_upper_triangular(u,na,T)
    % extract the upper triangular part of the visibility matrices for all time steps
    %
    % Args:
    %  u: (T,na,na)
    %  na: number of antennas
    %  T: number of time steps

    u2 = zeros(na*(na-1)/2,1);

    for t=1:T
        u_t = squeeze(u(t,:,:));
        u_tmp = u_t - triu(u_t);
        u_tmp = u_tmp(:);
        u_tmp = u_tmp(u_tmp~=0);
        u2((t-1)*na*(na-1)/2 + 1:t*na*(na-1)/2) = u_tmp;
    end

end