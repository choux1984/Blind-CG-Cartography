function [est_F,w,phi_col,evl_pnt] = estimate_F_and_w( s_check, Tx_pos , Rx_pos,  ini_F , myKfunc  , mu_w , mu_f, Nc , clustering_type, blind_ind , f,lambda_W)
%  
% 
% INPUT:
%   ini_F         Nx x Ny matrix with the initialization for the spatial loss field
%   sigma_2       Noise variance
%   t_slots       Number of measurements
%   myKfunc       Kernel function: must take a scalar ---> change to take
%                 two vectors
%   mu_w          Regularization parameter for weight function
%   mu_f          Regularization parameter for spatial loss field
%   Nc            Number of clusters
%                 Set Nc = [] for no clustering
%   blind_ind     1: Blind estimation
%                 0: Non-blind estimation
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
max_itr = 50;

% dependent variables
[N_x, N_y] = size(ini_F);
Ng = N_x * N_y; % # of grid points

t = length(s_check);
% Input check
assert( Nc <= Ng*t );

%%% estimation
[evl_pnt,idx_phi,phi_col]  = generate_feature_vecs(N_x,N_y,Tx_pos,Rx_pos,Nc,clustering_type,lambda_W);

% generation of K matrix %
K = myKmatrix(evl_pnt,myKfunc);

if blind_ind  % blind estimation
    [alpha,est_f] = blind_est(K,idx_phi,s_check,mu_w,ini_F(:),eps_error,max_itr,mu_f);
else  % non-blind
    alpha = est_alpha(K,idx_phi,ini_F(:),s_check,mu_w);
    est_f = est_field(K,idx_phi,alpha,s_check,Ng,R,mu_f,RK);
end

w = @(input)  myRepThm(alpha,evl_pnt,input,myKfunc);
est_F = reshape(est_f,N_x,N_y);

end