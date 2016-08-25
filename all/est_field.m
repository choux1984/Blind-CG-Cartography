function [est_f,A_K] = est_field(K,idx_phi,alpha,s_check,Ng,R,mu_f,RK)
    % Estimator for SLF
    
    % INPUT:
    %   K              Kernel matrix
    %   idx_phi        (n_measurements * Ng) x 1: Cluster indices where phi_i's belong to 
    %   alpha          (Ng-by-1) initialization of the spatial loss field (SLF)
    %   s_check        Noisy received signals
    %   Ng             Number of grid points
    %   t              Number of measurements
    %   R              N_c-by_N_c spatial covariance matrix (regularizer)
    %   mu             Tikhonov regularizer weight
    %
    % OUTPUT:
    %   est_f          Estimate of SLF
    %   A_K            t-by-Ng matrix
    t = length(s_check);
%     RK = compute_RK(K,idx_phi);
    temp_Ka = RK * alpha;
    A_K = reshape(temp_Ka,Ng,t)';
    est_f = (A_K' * A_K + t * mu_f * R)\ (A_K' * s_check);
end