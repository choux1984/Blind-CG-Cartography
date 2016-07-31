% Simulations for blind channel-gain cartography
% This is the main file
%
function simulations_BlindCG

mu_w_vec = 0.11;
K_std_vec = 0.24; %0.19
Nc_vec = 1000;
% lambda_vec = 0.2;
% K_std_vec = 0.15; %0.19
% Nc_vec = 2000;
t_slot_vec = 2400;

size_lambda = length(mu_w_vec);
size_K = length(K_std_vec);
size_Nc = length(Nc_vec);
size_t = length(t_slot_vec);
tic
for i = 1:size_lambda
    for j = 1 : size_K
        for k = 1 : size_Nc
            for t = 1 : size_t
                lambda = mu_w_vec(i);
                K_std = K_std_vec(j);
                Nc = Nc_vec(k);
                t_slots = t_slot_vec(t);
                simulation_01(lambda,K_std,Nc,t_slots);
            end
        end
    end
end

% simulation_02();
% simulation_03();
% simulation_04();
% simulation_05();
% simulation_06();

toc

end

% This simulation simply loads a SLF from a file and shows its estimate
% using the blind estimation algorithm
function simulation_01(mu_w,K_std,Nc,t_slots)

% Load SLF
close all
% filename1 = 'Map_6_6_v2.csv';  
% filename1 = 'Map_10_10_v2.csv';
filename1 = 'Map_30_30.csv';
F = csvread(filename1);
[N_x, N_y] = size(F);

% model parameters
sigma_2 = 0.0001;% noise variance
lambda_W = 3.5; % parameter to determine the threshold for nonzero weights
eps_W = 1.5 * 1e-1;    
w_fun = @(d,d_mu)  ( ((d./2).^2 + (d_mu).^2).^(-.5) )./ (pi.*d_mu);

% estimation parameters
mu_f = 1e-3;
myKfunc = @(input1,input2) exp(-norms(input1-input2).^2./(2 * K_std^2)); %myKfunc := kernel function
clustering_type = 'random';
blind_ind = 1;
if blind_ind
    ini_F = 20 * rand(N_x,N_y); 
%     ini_F = F; 

else
    ini_F = F;
end

% data generation
[s_check,Tx_pos,Rx_pos] = myRxSig(t_slots,F,sigma_2,w_fun,lambda_W,eps_W); %s_check := noisy received signal

% Estimation
[est_F,w_est,phi_col,evl_pnt] = estimate_F_and_w( s_check, Tx_pos , Rx_pos,  ini_F , myKfunc  , mu_w , mu_f,  Nc , clustering_type, blind_ind, F(:),lambda_W);

% est_f = (nW*nW' + t_slots * mu * eye(N_x*N_y))\(nW * s_check);
% est_F = reshape(est_f,N_x,N_y);

% Performance measurement
% MSE_F = norm(F-est_F)/numel(F)
% MSE_w = estimate_error_w(w_fun,w_est,phi_col)

% representation

% figure
% imagesc([F,est_F])
% colormap('gray');
% 
figure
imagesc(est_F)
colormap('gray');

a = imagesc(est_F);
imcontrast(a)

if est_F(10,6) < 0 
    est_F = -1 .* est_F;
    idx_sign = -1;
else
    idx_sign = 1;
end

h=figure
imagesc([est_F])
colormap('gray');
title(sprintf('Reconstructed field with K-std=%g, mu_w=%g, Nc=%g,t_slots=%g', K_std, mu_w,Nc,t_slots))
file_name=sprintf('est_f_K_std_%g_mu_w_%g_Nc_%g.fig',K_std, mu_w,Nc);
saveas(h,file_name)

compare_W(w_fun,w_est,N_x,N_y,lambda_W,eps_W,K_std,mu_w,Nc,idx_sign);

keyboard


end
