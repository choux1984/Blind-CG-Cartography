function [hat_alpha,hat_f] = blind_est(K,idx_phi,s_check,mu_w,ini_f,eps_error,max_itr,mu_f,resolution)
    % Estimators for alpha and F in a blind fashion
    
    % INPUT:
    %   K              Kernel matrix
    %   idx_phi        (n_measurements * Ng) x 1: Cluster indices where phi_i's belong to
    %   s_check        Noisy received signals
    %   mu_w           Parameter for Kernel smoothness regularizer
    %   t              Number of measurements = t_slots * num_CR_set
    %   ini_f          (Ng-by-1) initialization of the spatial loss field (SLF)
    %   eps_error      Tolerance of error btw current and previous iterates
    %   max_itr        Maximum iteration numbers
    %   mu_f           Regularization parameter for spatial loss field
    %   resolution     scaling factor to set a resolution of the reconstructed SLF (e.g. 1 = original resolution) 

    % OUTPUT:
    %   hat_alpha      Coeffcients of kernel functions parameterized by phi_i's 
    %                  [Non-clusetering] (t*Ng-by-t*Ng) matrix /
    %                  [Clustering] Nc-by-Nc
    %   hat_f          (Ng-by-1) estimated spatial loss field
    
    
    hat_f = ini_f;  
    
    Ng = length(ini_f);
    N_x = sqrt(Ng);
    N_y = N_x;
        
    TR = myiCovMat(20,20,resolution);
%     TR = eye(Ng); % Tikhonov reg.
    RK = compute_RK(K,idx_phi);
      
    % ADMM parameters
    rho = 1e-3;
    num_itr = 150;
    P = myPermMat(N_x,N_y);
    X = magic(N_x);
    D = TV_test(X);
    D1 = D(1:N_y*(N_x-1),:);
    
    for cnt = 1: max_itr
        %[S1]
        hat_alpha = est_alpha(K,idx_phi,hat_f,s_check,mu_w);
        %[S2]
%         [hat_f,A_K] = est_field(K,idx_phi,hat_alpha,s_check,Ng,TR,mu_f,RK);
%         [hat_f,A_K] = est_field_sp(K,idx_phi,hat_alpha,s_check,Ng,prev_hat_f,mu_f);
        prev_hat_f = hat_f;
        hat_f = est_field_ADMM(Ng,hat_alpha,mu_f,rho,num_itr,s_check,RK,D1,P,prev_hat_f);

        cnt
    end

end

function RK = compute_RK(K,idx_phi)
% Compute R*K. RK(i,:)= K(idx_phi(i),:) where R = Ng*t-by-Nc repetition
% matrix (cluster indicator)
%
% INPUT:
%   K             Nc-by-Nc Kernel matrix
%   idx_phi       Ng*t-by-1 Cluster indices that features belong to
%
% OUTPUT:
%   RK            R * K

RK = K(idx_phi,:);

end