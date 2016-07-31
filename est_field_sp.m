function [est_f,A_K] = est_field_sp(K,idx_phi,alpha,s_check,Ng,prev_hat_f,mu)
    % Estimator for SLF
    
    % INPUT:
    %   K              Kernel matrix
    %   idx_phi        (n_measurements * Ng) x 1: Cluster indices where phi_i's belong to 
    %   alpha          (Ng-by-1) initialization of the spatial loss field (SLF)
    %   s_check        Noisy received signals
    %   Ng             Number of grid points
    %   t              Number of measurements
    %   prev_hat_f     Previous estimate of f
    %   mu            Weight of l-1 norm regularizer 
    %
    % OUTPUT:
    %   est_f          Estimate of SLF
    %   A_K            t-by-Ng matrix
    t = length(s_check);
    Ng = size(prev_hat_f,1);
    
    RK = compute_RK(K,idx_phi);
    temp_Ka = RK * alpha;
    A_K = reshape(temp_Ka,Ng,t)';
    
    temp_tilde_s = s_check  - A_K * prev_hat_f;
    est_f = zeros(Ng,1);
    
    for l = 1 : Ng
        tilde_s = temp_tilde_s + A_K(:,l) * prev_hat_f(l);
        if norm(A_K(:,l)) == 0;
            est_f(l,1)=0;
        else
            func_var =  (tilde_s)' * A_K(:,l);
            est_f(l,1) = wthresh(func_var,'s', 0.5*t*mu)/(norm(A_K(:,l),2)^2);
        end
%         prev_hat_f(l) = est_f(l,1);
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