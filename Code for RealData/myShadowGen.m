function [s_check,omg_str1,Dist_str_dB] = myShadowGen(Pos_Tx,Pos_Rx,Pos_Tx_free,Pos_Rx_free,Ch_est_dB,Ch_est_dB_free,TxRxPair_file_name,TxRxPair_free_file_name)
% function to give an estimate of shadowing measurement from RSS

% Calculate the distances btw Rx's and Tx's
Dist_free = sqrt((Pos_Tx_free(:,1) - Pos_Rx_free(:,1)).^2 +(Pos_Tx_free(:,2) - Pos_Rx_free(:,2)).^2) * 0.3048;
Dist_str = sqrt((Pos_Tx(:,1) - Pos_Rx(:,1)).^2 +(Pos_Tx(:,2) - Pos_Rx(:,2)).^2) * 0.3048;
alpha = 2; % Compensation of pathloss : alpha = 2
Dist_free_dB = 10  * log(Dist_free)./log(10); % PL_free:=Pathloss of free space
Dist_str_dB = 10  * log(Dist_str)./log(10); % PL_str:=Pathloss of structure
PL_free = alpha * Dist_free_dB; % PL_free:=Pathloss of free space
PL_str = alpha * Dist_str_dB; % PL_str:=Pathloss of structure

% Estimation of Tx/Rx gain from free-space measurements
gain_free = Ch_est_dB_free + PL_free;
omg = csvread(TxRxPair_free_file_name);
omg_temp1 = omg(1:580,:);
omg_temp2 = omg(601:2400,:);
omg1 = [omg_temp1;omg_temp2];

gain = (omg1'*omg1 + 1e-12 * eye(size(omg1,2)))\(omg1'*gain_free); % Regularizer-based estimator

omg_str = csvread(TxRxPair_file_name);
omg_str_temp1 = omg_str(1:580,:);
omg_str_temp2 = omg_str(601:2400,:);
omg_str1 = [omg_str_temp1;omg_str_temp2];

gain_vec = omg_str1 * gain;
SH_comp_est = Ch_est_dB - gain_vec + PL_str;
temp_s1 = SH_comp_est(1:580,1);
temp_s2 = SH_comp_est(581:end,1);
s_check = [temp_s1;zeros(20,1);temp_s2];
% s_check = [temp_s1;temp_s2];


end