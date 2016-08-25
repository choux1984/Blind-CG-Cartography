%
%  FIGURES FOR THE PAPER ON BLIND CHANNEL GAIN CARTOGRAPHY
%
%  In this file, we essentially set parameters and invoke a simulator
%

classdef CGCartographySimulations < simFunctionSet
	
	properties
	
	end
	
	methods
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% %%  1. Sample simulations
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		% This is an empty simulation to check the environment
		function F = compute_fig_1001(obj,niter)
			
			T = 5
			
			F = [];
			
		end
			
		% Simple simulation where we draw a random figure
		function F = compute_fig_1002(obj,niter)
			
			x_values = 1:10;
			y_values =  randn(1,length(x_values));
			
			F = F_figure('X',x_values,'Y',y_values);
			
		end
			
		% Simple simulation where we draw a random figure and we label axes
		function F = compute_fig_1003(obj,niter)
			
			x_values = 1:10;
			y_values =  randn(3,length(x_values));
			
			F = F_figure('X',x_values,'Y',y_values,'leg',{'fun 1','fun 2','fun 3'},'xlab','X','ylab','FUNCTION VALUE');
			
		end
			
		% Simple simulation where we draw three random images
		function F = compute_fig_1004(obj,niter)
			
			F1 = F_figure('Z',rand(20,30));
			F2 = F_figure('Z',rand(20,30));
			F3 = F_figure('Z',rand(20,30));
			F = F_figure('multiplot_array',[F1 F2 F3]);

		end
			
		
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% %%  2. Simple non-blind simulations
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
		% This is a toy simulation
		function F = compute_fig_2001(obj,niter)
			
			
			% Create data genator
			m_F = csvread('Map_15_15.csv');
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			h_w = @(phi1,phi2) (1/sqrt(phi1))*(phi2<phi2+lambda_W/2); % normalized ellipse model
			s_measurementNum = 10;   
			s_noiseVar = 0.1;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar);
			
			% Create estimator
			ch_reg_f_type = 'tikhonov';
			mu_f = 1e-4;
			ini_F = randn(size(m_F));
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind');
			
			
			% SIMULATION
			% A) data generation
			[s_check,m_txPos,m_rxPos] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est] = est.estimate(s_check,m_txPos,m_rxPos)
			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%			
			
			% estimation parameters			
% 			mu_w = 0.26;
% 			
% 			K_std = 0.15; %0.19
% 			Nc = 10;
% 			rho = 1e-3;
% 			
% 			; % choose regularizer type for the SLF: 'T' = Tikhonov; 'S' = l_1; 'TV' = total variation						
% 			% myKfunc = @(input1,input2) exp(-norm(input1-input2).^2./(2 * K_std^2)); %myKfunc := kernel function
% 			% myKfunc = @(input1,input2) exp(-norm(input1-input2)./(2 * K_std^2)); %myKfunc := kernel function -> Expo. kernel			
% 			clustering_type = 'random';
			% initialization
			
			% data generation
			%[s_check,Tx_pos,Rx_pos] = myRxSig(t_slots,F,sigma_2,w_fun,lambda_W,eps_W); %s_check := noisy received signal
			
% 			% Estimation
 			%[est_F,w_est,phi_col,evl_pnt] = estimate_F_and_w( s_check, Tx_pos , Rx_pos,  ini_F , myKfunc  , mu_w , mu_f,  Nc , clustering_type, reg_f_type, blind_ind,lambda_W,rho);
% 			
% 			% representation
% 			
% 			
% 			h=figure
% 			imagesc([est_F])
% 			colormap('gray');
% 			
% 			save('est_F.mat','est_F')
% 			
% 			title(sprintf('Reconstructed field with K-std=%g, mu_f=%g, mu_w=%g, rho=%g, Nc=%g,t_slots=%g', K_std, mu_f, mu_w, rho, Nc,t_slots))
% 			file_name=sprintf('est_f_K_std_%g_mu_f_%g_mu_w_%g_rho_%g_Nc_%g.fig',K_std, mu_f, mu_w, rho, Nc);
% 			saveas(h,file_name)
% 			
% 			compare_W(w_fun,w_est,N_x,N_y,lambda_W,eps_W,K_std,mu_w,Nc,1,rho);
% 		
			
			F = [];
			
		end
		
		
        
	end
	
	
end
