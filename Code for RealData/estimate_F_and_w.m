function [est_F,w,phi_col,evl_pnt] = estimate_F_and_w( s_check, Tx_pos , Rx_pos,  ini_F , myKfunc  , mu_w , mu_f, Nc , clustering_type, blind_ind ,lambda_W,resolution)
%  
% 
% INPUT:
%   s_check       shadowing measurements
%   Tx_pos        t-by-2 Tx positions over a 2-D space
%   Rx_pos        t-by-2 Rx positions over a 2-D space
%   ini_F         Nx x Ny matrix with the initialization for the spatial loss field
%   myKfunc       Kernel function: must take a scalar ---> change to take
%                 two vectors
%   mu_w          Regularization parameter for weight function
%   mu_f          Regularization parameter for spatial loss field
%   Nc            Number of clusters
%                 Set Nc = [] for no clustering
%   blind_ind     1: Blind estimation
%                 0: Non-blind estimation
%   lambda_W      parameter to set a length of a semi-minor axis of an
%                 elliptical weight function
%   resolution    scaling factor to set a resolution of the reconstructed SLF (e.g. 1 = original resolution) 
%   Omega         (t-by-120) each row of Omega has two ones at the
%                 locations of Tx and Rx among 120 sensors
%   Dist_str_dB   distances between Tx/Rx in dB
%
% OUTPUT:
%   est_F         Nx x Ny matrix with the estimate of F
%   w             function of a 2D vector with the weight function
%                 if blind == 0, then w is the estimated weight function
%                 using the ground truth F
%   phi_col       2-by-(n_measurements*Ng): collection of all feature vectors

% constants
% num_CR_set = 1; % Number of secondary users in the CR network
eps_error = 1e-10;
max_itr = 200;

% dependent variables
[N_x, N_y] = size(ini_F);
Ng = N_x * N_y; % # of grid points

t = length(s_check);
% Input check
assert( Nc <= Ng*t );

%%% estimation
[evl_pnt,idx_phi,phi_col]  = generate_feature_vecs(N_x,N_y,Tx_pos,Rx_pos,Nc,clustering_type,lambda_W,resolution);

% generation of K matrix %
K = myKmatrix(evl_pnt,myKfunc);

if blind_ind  % blind estimation
    [alpha,est_f] = blind_est(K,idx_phi,s_check,mu_w,ini_F(:),eps_error,max_itr,mu_f,resolution);
else  % non-blind
    alpha = est_alpha(K,idx_phi,ini_F(:),s_check,mu_w);
    RK = K(idx_phi,:);
%     TR = myiCovMat(20,20,resolution);
%     [est_f,~] = est_field(K,idx_phi,alpha,s_check,Ng,TR,mu_f,RK);

%     est_f = est_field(K,idx_phi,alpha,s_check,Ng,R,mu_f,RK);
    
    % ADMM parameters
    rho = 1e-3; 
    num_itr = 200;
    P = myPermMat(N_x,N_y);
    dummy_F = magic(N_x);
    D = TV_test(dummy_F);
    D1 = D(1:N_y*(N_x-1),:); % for square F, D_x = D_y

    prev_hat_f = randn(Ng,1);
    est_f = est_field_ADMM(Ng,alpha,mu_f,rho,num_itr,s_check,RK,D1,P,prev_hat_f);

end

w = @(input)  myRepThm(alpha,evl_pnt,input,myKfunc);
est_F = reshape(est_f,N_x,N_y);

end