function simulation_01(mu_w,mu_f,K_std,Nc,t_slots,reg_f_type,filename,rho)

% Load SLF
close all
F = csvread(filename);
[N_x, N_y] = size(F);

% model parameters
sigma_2 = 0.1;% noise variance
lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
eps_W = 2 * 1e-2;    
w_fun = @(d,d_mu)  ( ((d./2).^2 + (d_mu).^2).^(-.5) )./ (pi.*d_mu);

% estimation parameters
% myKfunc = @(input1,input2) exp(-norm(input1-input2).^2./(2 * K_std^2)); %myKfunc := kernel function
myKfunc = @(input1,input2) exp(-norm(input1-input2)./(2 * K_std^2)); %myKfunc := kernel function -> Expo. kernel

clustering_type = 'random';
blind_ind = 1;
if blind_ind
    ini_F = 20 * rand(N_x,N_y); 
else
    ini_F = F;
end

% data generation
[s_check,Tx_pos,Rx_pos] = myRxSig(t_slots,F,sigma_2,w_fun,lambda_W,eps_W); %s_check := noisy received signal

% Estimation
[est_F,w_est,phi_col,evl_pnt] = estimate_F_and_w( s_check, Tx_pos , Rx_pos,  ini_F , myKfunc  , mu_w , mu_f,  Nc , clustering_type, reg_f_type, blind_ind,lambda_W,rho);

% representation
% figure
% imagesc([F,est_F])
% colormap('gray');

% a = imagesc(est_F);
% imcontrast(a)

% keyboard

h=figure
imagesc([est_F])
colormap('gray');

save('est_F.mat','est_F')

title(sprintf('Reconstructed field with K-std=%g, mu_f=%g, mu_w=%g, rho=%g, Nc=%g,t_slots=%g', K_std, mu_f, mu_w, rho, Nc,t_slots))
file_name=sprintf('est_f_K_std_%g_mu_f_%g_mu_w_%g_rho_%g_Nc_%g.fig',K_std, mu_f, mu_w, rho, Nc);
saveas(h,file_name)

compare_W(w_fun,w_est,N_x,N_y,lambda_W,eps_W,K_std,mu_w,Nc,1,rho);

% keyboard

end