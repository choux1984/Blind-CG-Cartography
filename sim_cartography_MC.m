function MSE = sim_cartography_MC( t_slots,F,sigma_2,w_fun , ini_F , myKfunc , lambda ,   Nc , clustering_type, blind_ind , n_iterations )
% Monte-Carlo Simulation for measuing MSE
%
% Input
%     t_slot             # of time slots
%     F                  Original spatial loss field
%     sigma_2            Noise variance
%     w_fun              Original weight function
%     ini_F              Initilization of a spatial loss field
%     myKfunc            Kernel function
%     lambda             Smoothness parameter
%     Nc                 Number of clusters
%     clustering_type    [kmeans] k-means; [random]: random sampling
%     blind_ind          [0]: Non-blind; [1]: Blind
%     n_iterations       Number of iterations for Monte-Carlo Simulation
%
% Output
%     MSE                Mean-squared error of ~~~

MSE_F_vec = zeros(1,n_iterations);
MSE_w_vec = zeros(1,n_iterations);

for n = 1:n_iterations

    % Data generation
    [s_check,Tx_pos,Rx_pos] = myRxSig(t_slots,F,sigma_2,w_fun)  ; %s_check := noisy received signal

    % Estimation
    [est_F,w_est,phi_col] = estimate_F_and_w(  s_check, Tx_pos , Rx_pos, ini_F , myKfunc , lambda ,   Nc , clustering_type, blind_ind); % step size

    % Measuruing performance 
    MSE_F_vec(n) = norm(F-est_F)/numel(F);
    MSE_w_vec(n) = estimate_error_w(w_fun,w_est,phi_col);
    
end

MSE = [mean(MSE_F_vec); mean(MSE_w_vec)];



end