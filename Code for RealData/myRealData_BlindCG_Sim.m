function myRealData_BlindCG_Sim
% INPUT
%  resolution         integer to decide the resol. of the discretized SLF
close all
resolution= 1.25;

% Load files
filename1 = 'Tx_position.csv'; % v2 for non-missing case (step6)
filename2 = 'Rx_position.csv'; 
filename3 = 'total_data.csv';
filename4 = 'Tx_position_free.csv'; filename5 = 'Rx_position_free.csv';
filename6 = 'total_data_free.csv';
filename7 = 'omg.csv'; filename8 = 'omg1.csv';
filename9 = 'est_F_0.csv';

Tx_Pos = csvread(filename1); Rx_Pos = csvread(filename2); Ch_est = csvread(filename3); % Structured scenario
Tx_Pos_free = csvread(filename4);Rx_Pos_free = csvread(filename5);% free-space scenario
Ch_est_free_temp= csvread(filename6);Ch_est_free = [Ch_est_free_temp(1:580,1);Ch_est_free_temp(601:end,1)]; % Channel gains
% Convert Channel_est to dB scaled values
Ch_est_dB = 10 * log(Ch_est) / log(10);
Ch_est_dB_free = 10 * log(Ch_est_free) / log(10);

[Tx_Pos_total,Rx_Pos_total,Tx_Pos_total_free,Rx_Pos_total_free,temp_Tx_Pos_total,temp_Rx_Pos_total] = myTotalPosition(Tx_Pos,Rx_Pos,Tx_Pos_free,Rx_Pos_free);

% Estimate shadowing measurements
[s_check,omg_str,Dist_str_dB] = myShadowGen(Tx_Pos_total,Rx_Pos_total,Tx_Pos_total_free,Rx_Pos_total_free,Ch_est_dB,Ch_est_dB_free,filename7,filename8);

% estimation parameters
N_x = 20 * resolution + 1; % Size of the discretized SLF (N_x-by-N_y)
N_y = 20 * resolution + 1;
mu_f = 0.05;
mu_w = 1e-2;
Nc = 5;
K_std = 0.05;
myKfunc = @(input1,input2) exp(-norms(input1-input2).^2./(2 * K_std^2)); %myKfunc := kernel function
clustering_type = 'random';
blind_ind = 0;
% ini_F = randn(N_x,N_y); 
ini_F = csvread(filename9);
lambda_W = 0.1224; % parameter to determine the threshold for nonzero weights
w_fun = @(d,d_mu)  ( ((d./2).^2 + (d_mu).^2).^(-.5) )./ (pi.*d_mu);
eps_W = 1e-2;    

% Estimation
[est_F,w_est,phi_col,evl_pnt] = estimate_F_and_w( s_check, temp_Tx_Pos_total', temp_Rx_Pos_total',  ini_F , myKfunc  , mu_w , mu_f,  Nc , clustering_type, blind_ind,lambda_W,resolution);

figure
imagesc(est_F)
colormap('gray');

a = imagesc(est_F);
imcontrast(a)

if est_F(6,8) < 0 
    est_F = -1 .* est_F;
    idx_sign = -1;
else
    idx_sign = 1;
end

% h=figure
% imagesc([est_F])
% colormap('gray');
% title(sprintf('Reconstructed field with K-std=%g, mu_w=%g, Nc=%g,t_slots=%g', K_std, mu_w,Nc,t_slots))
% file_name=sprintf('est_f_K_std_%g_mu_w_%g_Nc_%g.fig',K_std, mu_w,Nc);
% saveas(h,file_name)
% keyboard
compare_W(w_fun,w_est,N_x,N_y,lambda_W,eps_W,K_std,mu_w,Nc,idx_sign);

keyboard

end

function [Tx_Pos_total,Rx_Pos_total,Tx_Pos_total_free,Rx_Pos_total_free,temp_Tx_Pos_total,temp_Rx_Pos_total] = myTotalPosition(Tx_Pos,Rx_Pos,Tx_Pos_free,Rx_Pos_free)
% Preprocess to aggregate sensor positions according to the access order of Rx00-09 from Tx00-09 % 

for i = 1 : 24    
    temp_Tx = Tx_Pos(1 + (i-1)*10: i*10,:);
    temp_Tx_free = Tx_Pos_free(1 + (i-1)*10: i*10,:);
    temp_Rx = Rx_Pos(1 + (i-1)*10: i*10,:);
    temp_Rx_free = Rx_Pos_free(1 + (i-1)*10: i*10,:);
    for j = 1 : 10     
        temp_Tx_Pos_total(1 + (i-1)*100 + (j-1)*10:(i-1)*100 + j*10,:) = temp_Tx;
        temp_Rx_Pos_total(1 + (i-1)*100 + (j-1)*10:(i-1)*100 + j*10,:) = ones(10,1) * temp_Rx(j,:);
        temp_Tx_Pos_total_free(1 + (i-1)*100 + (j-1)*10:(i-1)*100 + j*10,:) = temp_Tx_free;
        temp_Rx_Pos_total_free(1 + (i-1)*100 + (j-1)*10:(i-1)*100 + j*10,:) = ones(10,1) * temp_Rx_free(j,:);
    end
end

% Remove the tx/rx positions from 580-600 (unobserved due to error)
Tx_Pos_total = [temp_Tx_Pos_total(1:580,:);temp_Tx_Pos_total(601:end,:)];
Rx_Pos_total = [temp_Rx_Pos_total(1:580,:);temp_Rx_Pos_total(601:end,:)];
Tx_Pos_total_free = [temp_Tx_Pos_total_free(1:580,:);temp_Tx_Pos_total_free(601:end,:)];
Rx_Pos_total_free = [temp_Rx_Pos_total_free(1:580,:);temp_Rx_Pos_total_free(601:end,:)];

end