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
		% %%  2. Simple non-blind simulations with synthetic data
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
		% This is a toy simulation for a non-blind shadow loss field
		% estimation. No simultaneous calibration (estimator knows path
		% loss exponent and sensor gains)
		function F = compute_fig_2001(obj,niter)
			
			
			% Create data genator
			m_F = csvread('Map_15_15.csv');
			m_F = m_F/max(max(m_F));
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			h_w = @(phi1,phi2) (1/sqrt(phi1))*(phi2<phi1+lambda_W/2); % normalized ellipse model
			s_measurementNum = 200;   
			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar);
			
			% Create estimator
			ch_reg_f_type = 'tikhonov'; %'l1_PCO'; %'totalvariation';
			mu_f = 1e-4;
			ini_F = randn(size(m_F));
            rho =  1e-2;
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','rho',rho);
			
			
			% SIMULATION
			% A) data generation
			[s_check,m_txPos,m_rxPos] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est] = est.estimate(s_check,m_txPos,m_rxPos);
			
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');			
			F = F_figure('multiplot_array',[F1 F2]);
			
			
		
			
		end% This is a toy simulation for a non-blind shadow loss field
		
		
		% comparison of different regularizers (to be coded)
		
		
		
		% non-blind with simultaneous calibration
		
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% %%  3. Simple blind simulations with synthetic data
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
		% Simulation with blind estimation
		function F = compute_fig_3001(obj,niter)
			
			F = [];
		end
		
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% %%  4. Simple non-blind simulations with REAL data
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		% without simultaneous calibration
		
		
		
		
		% with simultaneous calibration
		
				
        
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% %%  5. Simple blind simulations with REAL data
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
		
	end
	
	
end
