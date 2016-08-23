function est_f = est_field_ADMM(Ng,alpha,mu_f,rho,num_itr,s_check,RK,D,P,f)
    % Estimator for SLF via ADMM
    
    % INPUT:
    %   Ng             Number of grid points
    %   alpha          (Ng-by-1) initialization of the spatial loss field (SLF)
    %   mu_f           TV regularizer weight
    %   rho            step size for dual ascents
    %   num_itr        Number of iterations for ADMM
    %   s_check        Noisy received signals
    %     D            N_y(N_x-1)-by-Ng discrete gradient matrix
    %     P            Permutation matrix to make row-major vector of 'f'  
    %     f            Initialization of f
    %
    % OUTPUT:
    %   est_f          Estimate of SLF

    T = length(s_check);
    temp_Ka = RK * alpha;
    A = reshape(temp_Ka,Ng,T)';    

    [N_Dx,N_g] = size(D);
    
    til_A = (2/T) * A' * A;
    til_D = D'*D;
    DP = D * P;
    til_P = DP' * DP;
    til_H = (til_A + rho * (til_D + til_P) );
    til_s = (2/T) * A' * s_check;
    
    % Initialize
%     f = randn(N_g,1);
    gamma_x = zeros(N_Dx,1);
    gamma_y = zeros(N_Dx,1);
    d_x = zeros(N_Dx,1);
    d_y = zeros(N_Dx,1);
    
    for k = 1 : num_itr
%         k
        % [S1] Dual ascent for aux variables
        gamma_x = gamma_x + rho * (D*f - d_x);
        gamma_y = gamma_y + rho * (DP*f - d_y);  
        
        % [S2] Updates of aux variables d_x and d_y
        d_x = wthresh(D*f + gamma_x./rho,'s', mu_f/rho);
        d_y = wthresh(DP*f + gamma_y./rho,'s', mu_f/rho);
        
        % [S3] Updates of SLF        
        f = ( til_H ) \(D' *(rho * d_x - gamma_x) + DP' *(rho * d_y -gamma_y) + til_s);
%         k
%         cost(k) = (1/T) * norm(A*f - s_check,2)^2 + mu_f *(sum(abs(D*f)) + sum(abs(DP*f)));
    end
    est_f = f;
%     figure
%     plot(cost);
    
end