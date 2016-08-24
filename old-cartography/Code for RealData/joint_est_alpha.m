function [TR_gain,gamma,alpha] = joint_est_alpha(K,idx_phi,f,Ch_est_dB,mu_w,Omega,Dist_str_dB,til_K)
    % Estimator for alpha, tx/rx gains, and path-loss exponent
    % the number of sensors = 120
    
    % INPUT:
    %   K              Kernel matrix
    %   idx_phi        (n_measurements * Ng) x 1: Cluster indices where phi_i's belong to 
    %                  idx_phi(1) means that the location of 1 on the 1st row of R matrix 
    %   f              (Ng-by-1) initialization of the spatial loss field (SLF)
    %   Ch_est_dB      (t-by-1) CG measurements in dB
    %   mu_w           Parameter for Kernel smoothness regularizer
    %   Omega          (t-by-120) each row of Omega has two ones at the
    %                  locations of Tx and Rx among 120 sensors
    %   Dist_str_dB    distances between Tx/Rx in dB
    %   til_K          (120 + Nc + 1) Block diagonal matrix where the last block is 'K'
    %                  and the others are zero matrices
    %
    % OUTPUT:
    %   alpha          Estimate of alpha
    
    t = size(Ch_est_dB,1);
    Nc = size(K,2);

    temp_fR = myfRMat(t,Nc,f,idx_phi);    
    K1 = temp_fR * K;
    
    til_Omega = [Omega,-1 * Dist_str_dB, -1 * K1];
    
    A = (til_Omega' * til_Omega + mu_w * t * til_K);
    til_alpha = A\(til_Omega' * Ch_est_dB);
    
    TR_gain = til_alpha(1:120,1);
    gamma = til_alpha(121,1);
    alpha = til_alpha(122:end,1);
    
%     A = (temp_fR' * K1 + mu_w * t * eye(Nc));
%     alpha = A\(temp_fR' * s_check);    
end

function fR_Mat = myfRMat(t,Nc,f,idx_phi)
% fR_mat = ( I \kron f^T ) R
% i) R_tilde:=[R_1' ; ... ; R_T'] * f
% ii) Rf_vec := sum(R_tilde,2)
% iii) fR_Mat = unvec(Rf_vec)'

Ng = size(f,1);

idx_phi_mat = reshape(idx_phi,Ng,t);
row_inds = idx_phi_mat + Nc*ones(Ng,1)*(0:t-1);
col_inds = repmat(1:Ng,1,t)';
entry_vals = repmat(f,t,1);

Rf_vec = full(sum(sparse(row_inds(:),col_inds,entry_vals,Nc*t,Ng,length(idx_phi)),2));

fR_Mat = reshape(Rf_vec,Nc,t)';


end