% INPUT
%  resolution         integer to decide the resol. of the discretized SLF!!
resolution= 2;

% Load files
filename1 = 'Tx_position.csv'; % v2 for non-missing case (step6)
filename2 = 'Rx_position.csv'; 
filename3 = 'total_data.csv';
filename4 = 'Tx_position_free.csv'; filename5 = 'Rx_position_free.csv';
filename6 = 'total_data_free.csv';
filename7 = 'omg.csv'; filename8 = 'omg1.csv';
filename9 = 'est_F_0_42_by_42.csv';

Tx_Pos = csvread(filename1); Rx_Pos = csvread(filename2); Ch_est = csvread(filename3); % Structured scenario
Tx_Pos_free = csvread(filename4);Rx_Pos_free = csvread(filename5);% free-space scenario
Ch_est_free_temp= csvread(filename6);Ch_est_free = [Ch_est_free_temp(1:580,1);Ch_est_free_temp(601:end,1)]; % Channel gains
% Convert Channel_est to dB scaled values
Ch_est_dB = 10 * log(Ch_est) / log(10);
Ch_est_dB_free = 10 * log(Ch_est_free) / log(10);

[Tx_Pos_total,Rx_Pos_total,Tx_Pos_total_free,Rx_Pos_total_free,temp_Tx_Pos_total,temp_Rx_Pos_total] = myTotalPosition(Tx_Pos,Rx_Pos,Tx_Pos_free,Rx_Pos_free);

% Estimate shadowing measurements
[s_check,Omega,Dist_str_dB]  = myShadowGen(Tx_Pos_total,Rx_Pos_total,Tx_Pos_total_free,Rx_Pos_total_free,Ch_est_dB,Ch_est_dB_free,filename7,filename8);

% estimation parameters
N_x = 20 * resolution + 1; % Size of the discretized SLF (N_x-by-N_y)
N_y = 20 * resolution + 1;
mu_f = 1e-2; % 1e-1
mu_w = 2 * 1e-3; % 2 * 1e-3
Nc = 3000;
K_std = 0.2;
% myKfunc = @(input1,input2) exp(-norm(input1-input2).^2./(2 * K_std^2)); %myKfunc := kernel function
myKfunc = @(input1,input2) exp(-norm(input1-input2)./(2 * K_std^2)); %myKfunc := kernel function

clustering_type = 'random';
blind_ind = 1;
ini_F = 2 * rand(N_x,N_y) - 1; 
% ini_F = csvread(filename9);
lambda_W = 0.1224; % parameter to determine the threshold for nonzero weights
w_fun = @(d,d_mu)  ( ((d./2).^2 + (d_mu).^2).^(-.5) )./ (pi.*d_mu);
eps_W = 1e-3;    

% Estimation
% [est_F,w_est,phi_col,evl_pnt] = est_joint_F_w( Ch_est_dB, Tx_Pos_total', Rx_Pos_total',  ini_F , myKfunc  , mu_w , mu_f, Nc , clustering_type, blind_ind ,lambda_W,resolution,Omega,Dist_str_dB);
[est_F,w_est,phi_col,evl_pnt] = estimate_F_and_w( s_check, temp_Tx_Pos_total', temp_Rx_Pos_total',  ini_F , myKfunc  , mu_w , mu_f,  Nc , clustering_type, blind_ind,lambda_W,resolution);

h=figure
imagesc(est_F)
colormap('gray');

save('est_F.mat','est_F')

a = imagesc(est_F);
imcontrast(a)

title(sprintf('Reconstructed field with K-std=%g, mu_f=%g, mu_w=%g, Nc=%g,t_slots=%g', K_std, mu_f,mu_w,Nc))
file_name=sprintf('est_f_K_std_%g_mu_f_%g_mu_w_%g_Nc_%g.fig',K_std, mu_f, mu_w,Nc);
saveas(h,file_name)

compare_W(w_fun,w_est,N_x,N_y,lambda_W,eps_W,K_std,mu_w,Nc);

