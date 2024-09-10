function [DG, adjoint_DG] = util_precompute_op(G, param_ROP)
    % Compte the dense matrix composing the compression operator D with the off-grid interpolation operator G. 
    % G : interpolation matrix (eventually with the visibility weighting W)
    % param_ROP : struct with the parameters for the rank-one projections.

    % Options are:
    % (1) sparse_matrix
    % (2) batchwise
    % (3) for_loop
    precomp_type = 1;

    if precomp_type == 1
        % (1) Precompute from a product of sparse matrices

        alpha = param_ROP.alpha;
        beta = param_ROP.beta;
        Gamma = param_ROP.Gamma;

        D = sep_ROP_sparse_mat(alpha, beta);

        whos G;
        whos D;

        tic;
        DG = D*G;
        toc;
        whos DG

        clear D;
        
        M = modul_ROP_sparse_mat(Gamma, param_ROP.Np);

        whos M;

        tic; 
        MDG = M * DG;
        toc;
        whos MDG

        clear M;

    elseif precomp_type == 3

        % (2) Precompute by iterating over the batches

    elseif precomp_type == 3
        % (3) Precompute by looping over the columns of G

        % % Make histogram of the columns of G
        % figure;
        % histogram(sum(abs(G),1));
        % title('Histogram of the columns of G');
        % hold on;

        % Select indices of the nonzero columns in G with a tolerance
        tol = 1e-6;
        ind = find(sum(abs(G),1)>tol);
        % size(G)
        % length(ind)

        % Preallocate the matrix
        Nmeas = param_ROP.Np * param_ROP.Nm;
        MDG = zeros(Nmeas,length(ind));

        % Fill the matrix
        tic;
        for i=1:length(ind)
            if mod(i,100)==0
                disp(['Processing column ',num2str(i),' of ',num2str(length(ind))]);
            end
            MDG(:,i) = preop(G(:,ind(i)));
        end
        toc;

        tic;
        ncol = 5;
        count=0;
        while count*ncol < length(ind)
            disp(['Processing columns: ',num2str(count*ncol/length(ind)*100),'%']);
            MDG(:,count*ncol+1:(count+1)*ncol) = preop(G(:,count*ncol+1:(count+1)*ncol));
            count = count + 1;
        end
        MDG = preop(G(:,ind));
        MDG = MDG(:,ind);
        toc;

    end


    save('MDG.mat','MDG');

    figure;
    imshow(abs(MDG),[]);
    title('Matrix of the precomputed ROP operator');
    colorbar;

end