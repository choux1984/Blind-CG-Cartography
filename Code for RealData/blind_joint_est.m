function [hat_alpha,hat_f] = blind_joint_est(K,idx_phi,Ch_est_dB,mu_w,mu_f,ini_f,max_itr,resolution,Omega,Dist_str_dB)
    % Estimators for alpha and F in a blind fashion
    
    % INPUT:
    %   K              Kernel matrix
    %   idx_phi        (n_measurements * Ng) x 1: Cluster indices where phi_i's belong to
    %   Ch_est_dB      CG measurements in dB
    %   mu_w           Parameter for Kernel smoothness regularizer
    %   mu_w           Parameter for SLF regularizer
    %   ini_f          (Ng-by-1) initialization of the spatial loss field (SLF)
    %   max_itr        Maximum iteration numbers
    %   resolution     scaling factor to set a resolution of the reconstructed SLF (e.g. 1 = original resolution) 
    %   Omega          (t-by-120) each row of Omega has two ones at the
    %                  locations of Tx and Rx among 120 sensors
    %   Dist_str_dB    distances between Tx/Rx in dB
    % OUTPUT:
    %   hat_alpha      Coeffcients of kernel functions parameterized by phi_i's 
    %                  [Non-clusetering] (t*Ng-by-t*Ng) matrix /
    %                  [Clustering] Nc-by-Nc
    %   hat_f          (Ng-by-1) estimated spatial loss field
    
    
    hat_f = ini_f;  
    
    Ng = length(ini_f);
    N_x = sqrt(Ng);
    N_y = N_x;
        
%     TR = myiCovMat(20,20,resolution);

    TR = eye(Ng); % Tikhonov reg.
    RK = compute_RK(K,idx_phi);
    
    dim_g = 120 + 1 + size(K,2);
    til_K = zeros(dim_g,dim_g);
    til_K(1:121,1:121) = 1e-10 * eye(121);
    til_K(122:end,122:end) = K;
    
    % ADMM parameters
    rho = 1e-4;
    num_itr = 50;
    P = myPermMat(N_x,N_y);
    X = magic(N_x);
    D = TV_test(X);
    D1 = D(1:N_y*(N_x-1),:);
    
    for cnt = 1: max_itr
        %[S1]
%         hat_alpha = est_alpha(K,idx_phi,hat_f,s_check,mu_w);
        [TR_gain,gamma,hat_alpha] = joint_est_alpha(K,idx_phi,hat_f,Ch_est_dB,mu_w,Omega,Dist_str_dB,til_K);
        s_check = Omega * TR_gain - Dist_str_dB * gamma - Ch_est_dB;
        %[S2]
        [hat_f,A_K] = est_field(K,idx_phi,hat_alpha,s_check,Ng,TR,mu_f,RK);
%         [hat_f,A_K] = est_field_sp(K,idx_phi,hat_alpha,s_check,Ng,prev_hat_f,mu_f);
%         prev_hat_f = hat_f;
%         hat_f = est_field_ADMM(Ng,hat_alpha,mu_f,rho,num_itr,s_check,RK,D1,P,prev_hat_f);

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