function [est_F,w,phi_col,evl_pnt] = est_joint_F_w( Ch_est_dB, Tx_pos , Rx_pos,  ini_F , myKfunc  , mu_w , mu_f, Nc , clustering_type, blind_ind ,lambda_W,resolution,Omega,Dist_str_dB)
%  
% 
% INPUT:
%   Ch_est_dB        CG measurements in dB
%   Tx_pos           t-by-2 Tx positions over a 2-D space
%   Rx_pos           t-by-2 Rx positions over a 2-D space
%   ini_F            Nx x Ny matrix with the initialization for the spatial loss field
%   myKfunc          Kernel function: must take a scalar ---> change to take
%                    two vectors
%   mu_w             Regularization parameter for weight function
%   mu_f             Regularization parameter for spatial loss field
%   Nc               Number of clusters
%                    Set Nc = [] for no clustering
%   clustering_type 'random': random sampling
%                   'kmeans': kmeans  
%   blind_ind        1: Blind estimation
%                    0: Non-blind estimation
%   lambda_W         parameter to set a length of a semi-minor axis of an
%                    elliptical weight function
%   resolution       scaling factor to set a resolution of the reconstructed SLF (e.g. 1 = original resolution) 
%   Omega            (t-by-120) each row of Omega has two ones at the
%                    locations of Tx and Rx among 120 sensors
%   Dist_str_dB      distances between Tx/Rx in dB
%
% OUTPUT:
%   est_F         Nx x Ny matrix with the estimate of F
%   w             function of a 2D vector with the weight function
%                 if blind == 0, then w is the estimated weight function
%                 using the ground truth F
%   phi_col       2-by-(n_measurements*Ng): collection of all feature vectors
%   evl_pnt       2-by-Nc mat: collections of testing points (e.g. centroids)


% constants
max_itr = 1;

% dependent variables
[N_x, N_y] = size(ini_F);
Ng = N_x * N_y; % # of grid points
t = length(Ch_est_dB);

% Input check
assert( Nc <= Ng*t );

%%% estimation
[evl_pnt,idx_phi,phi_col]  = generate_feature_vecs(N_x,N_y,Tx_pos,Rx_pos,Nc,clustering_type,lambda_W,resolution);

% generation of K matrix %
K = myKmatrix(evl_pnt,myKfunc);

if blind_ind  % blind estimation
    [alpha,est_f] = blind_joint_est(K,idx_phi,Ch_est_dB,mu_w,mu_f,ini_F(:),max_itr,resolution,Omega,Dist_str_dB);
else  % non-blind    
    RK = K(idx_phi,:);
    dim_g = 120 + 1 + size(K,2);
    til_K = zeros(dim_g,dim_g);
    til_K(1:121,1:121) = 1e-10 * eye(121);
    til_K(122:end,122:end) = K;
    
    [TR_gain,gamma,alpha] = joint_est_alpha(K,idx_phi,ini_F(:),Ch_est_dB,mu_w,Omega,Dist_str_dB,til_K);
    s_check = Omega * TR_gain - Dist_str_dB * gamma - Ch_est_dB;
    
    TR = myiCovMat(20,20,resolution);
    [est_f,~] = est_field(K,idx_phi,alpha,s_check,Ng,TR,mu_f,RK);

%     est_f = est_field(K,idx_phi,alpha,s_check,Ng,R,mu_f,RK);
    
%     % ADMM parameters
%     rho = 1e-4; 
%     num_itr = 200;
%     P = myPermMat(N_x,N_y);
%     dummy_F = magic(N_x);
%     D = TV_test(dummy_F);
%     D1 = D(1:N_y*(N_x-1),:); % for square F, D_x = D_y
% 
%     prev_hat_f = randn(Ng,1);
%     est_f = est_field_ADMM(Ng,alpha,mu_f,rho,num_itr,s_check,RK,D1,P,prev_hat_f);

end

w = @(input)  myRepThm(alpha,evl_pnt,input,myKfunc);
est_F = reshape(est_f,N_x,N_y);

end