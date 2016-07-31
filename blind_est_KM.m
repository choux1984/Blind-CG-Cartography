function [hat_alpha,hat_f] = blind_est_KM(K,idx_phi,s_check,lambda,ini_f,eps_error,max_itr,gamma)
    % Estimators for alpha and F in a blind fashion
    
    % INPUT:
    %   K              Kernel matrix
    %   s_check        Noisy received signals
    %   lambda         Parameter for Kernel smoothness regularizer
    %   t              Number of measurements = t_slots * num_CR_set
    %   ini_f          (Ng-by-1) initialization of the spatial loss field (SLF)
    %   eps_error      Tolerance of error btw current and previous iterates
    %   max_itr        Maximum iteration numbers
    %   gamma          
    % OUTPUT:
    %   hat_alpha      Coeffcients of kernel functions parameterized by phi_i's 
    %                  [Non-clusetering] (t*Ng-by-t*Ng) matrix /
    %                  [Clustering] Nc-by-Nc
    %   hat_f          (Ng-by-1) estimated spatial loss field
    
    

    hat_f = ini_f;
%     hat_alpha = rand(size(K,1),1);
%     gamma_alp = 1;

    Ng = length(ini_f);
    t = length(s_check);
    
    debug = 0;         
    if debug
        prev_est_f = inf;
        prev_est_alpha = inf;
    end

    for cnt = 1: max_itr
               
        %[S1]
%         prev_hat_alpha = hat_alpha;
        hat_alpha = est_alpha(K,idx_phi,hat_f,s_check,lambda);
%         hat_alpha = gamma_alp * hat_alpha + (1 - gamma_alp) * prev_hat_alpha;       
        %[S2]
        prev_hat_f = hat_f;
        [hat_f,A_K] = est_field(K,idx_phi,hat_alpha,s_check,Ng);
        hat_f = gamma * hat_f + (1 - gamma) * prev_hat_f;
        
        
    end
    

end