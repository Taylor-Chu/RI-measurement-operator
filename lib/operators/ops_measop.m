function [measop, adjoint_measop, y] = ops_measop(vis, G, Ft, IFt, param_weighting, tau, param_ROP)
    % Generate the measurement op and its adjoint from a sampling pattern and
    % user input settings
    % operator (adapted from original code associated with
    % :cite:p:`Fessler2003`).
    %
    % Parameters
    % ----------
    %
    % G : function handle
    %     Visibility interpolation operator
    % Ft : function handle
    %     Fourier transform operator
    % IFt : function handle
    %     Inverse Fourier transform operator
    % param_weighting : struct
    %     weighting_on : bool
    %         Structure containing the parameters for visibility weighting
    %     nWimag : double[:]
    %         weights inferred from the density of the sampling (uniform/Briggs)
    % tau : float
    %     Noise level
    % param_ROP : struct
    %     Structure containing the parameters for applying ROPs on the measurements
    %
    % Returns
    % -------
    % measop : function handle
    %     Measurement operator
    % adjoint_measop : function handle
    %     Adjoint of the measurement operator
    % y : double[:]
    %     Measurement vector
    
    nWimag = param_weighting.nWimag;
    weighting_on = param_weighting.weighting_on;

    % (optionally) apply visibility weighting
    if weighting_on
        % nW = (1 / tau) * ones(na^2*nTimeSamples,1);
        nW = (1 / tau) * nWimag;
        [W, Wt] = op_vis_weighting(nW);
        vis = W(vis);
    end
    
    % (optionally) apply ROPs
    if param_ROP.use_ROP

        %% compute ROP operator
        [D, Dt] = op_ROP(param_ROP);

        y = D(vis);

        % Precompute forward operator
        if strcmp(param_ROP.type, 'modul')

            % [D, ~, M] = op_ROP(param_ROP);

            % if weighting_on
            %     preop = @(x) ( D(W(x)) );
            % else 
            %     preop = @(x) ( D(x) );
            % end

            if weighting_on
                measop = @(x) ( D(W(G * Ft(x))) ) ; 
                adjoint_measop = @(y) real(IFt(G' * Wt(Dt(y))));
            else 
                measop = @(x) ( D(G * Ft(x)) ) ; 
                adjoint_measop = @(y) real(IFt(G' * Dt(y)));
            end
            
            % Precompute as a matrix

            % % Make histogram of the columns of G
            % size(G)
            % figure;
            % size(sum(abs(G),1))
            % histogram(sum(abs(G),1));
            % title('Histogram of the columns of G');
            % hold on;

            % % Select indices of the nonzero columns in G with a tolerance
            % tol = 1e-6;
            % ind = find(sum(abs(G),1)>tol);
            % length(ind)

            % % Preallocate the matrix
            % Nmeas = param_ROP.Np * param_ROP.Nm
            % MDG = zeros(Nmeas,length(ind));

            % % Fill the matrix
            % tic;
            % for i=1:length(ind)
            %     if mod(i,100)==0
            %         disp(['Processing column ',num2str(i),' of ',num2str(length(ind))]);
            %     end
            %     MDG(:,i) = preop(G(:,ind(i)));
            % end
            % toc;

            % tic;
            % ncol = 5;
            % count=0;
            % while count*ncol < length(ind)
            %     disp(['Processing columns: ',num2str(count*ncol/length(ind)*100),'%']);
            %     MDG(:,count*ncol+1:(count+1)*ncol) = preop(G(:,count*ncol+1:(count+1)*ncol));
            %     count = count + 1;
            % end
            % MDG = preop(G(:,ind));
            % MDG = MDG(:,ind);
            % toc;

            % tic;
            % DG = D*G;
            % toc;
            % whos DG

            % tic; 
            % MDG = M * DG;
            % toc;
            % whos MDG

            % save('MDG.mat','MDG');

            % figure;
            % imshow(abs(MDG),[]);
            % title('Matrix of the precomputed ROP operator');
            % colorbar;

            % measop = @(x) ( MDG * Ft(x) ) ;
            % adjoint_measop = @(y) real(IFt(MDG' * y));
            
        else
            % %% compute ROP operator
            % [D, Dt] = op_ROP(param_ROP);

            if weighting_on
                measop = @(x) ( D(W(G * Ft(x))) ) ; 
                adjoint_measop = @(y) real(IFt(G' * Wt(Dt(y))));
            else 
                measop = @(x) ( D(G * Ft(x)) ) ; 
                adjoint_measop = @(y) real(IFt(G' * Dt(y)));
            end
        end

    else
        if weighting_on
            measop = @(x) ( W(G * Ft(x)) ) ; 
            adjoint_measop = @(y) real(IFt(G' * Wt(y)));
        else 
            measop = @(x) ( G * Ft(x) ) ; 
            adjoint_measop = @(y) real(IFt(G' * y));
        end

        y = vis;
    end
    
end    