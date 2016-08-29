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
		% loss exponent and sensor gains).  
		function F = compute_fig_2001(obj,niter)
						
			% Create data genator
			m_F = csvread('Map_15_15.csv');
			m_F = m_F/max(max(m_F));
			s_numGrid = size(m_F,1) * size(m_F,2);
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			h_w = @(phi1,phi2) (1/sqrt(phi1))*(phi2<phi1+lambda_W/2); % normalized ellipse model
			s_measurementNum = 200;
            s_pathLossExponent = 2;
            v_gains = [];
			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
			
			% Create estimator
			ch_reg_f_type = 'tikhonov'; %'l1_PCO'; %'totalvariation';
			m_tikhonov = eye(s_numGrid);
			ch_calibrationType = 'none'; % 'none','previous','simultaneous'
			mu_f = 1e-4;
			ini_F = randn(size(m_F));
            rho =  1e-2; % for ISTA, rho is roughtly 2.5
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'m_tikhonov',m_tikhonov);
						
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements,[]);
			
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');			
			F = F_figure('multiplot_array',[F1 F2]);		
			
        end
		
		
		% Comparison of different regularizers. WLOG, assume that an estimator knows path
		% loss exponent and sensor gains.  
		function F = compute_fig_2002(obj,niter)
						
			% Create data genator
			m_F = csvread('Map_15_15.csv');
			m_F = m_F/max(max(m_F));
			s_numGrid = size(m_F,1) * size(m_F,2);
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			h_w = @(phi1,phi2) (1/sqrt(phi1))*(phi2<phi1+lambda_W/2); % normalized ellipse model
			s_measurementNum = 200;
            s_pathLossExponent = 2;
            v_gains = [];
			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
			
			% Create estimator
			m_tikhonov = eye(s_numGrid);
			ch_calibrationType = 'none'; % 'none','previous','simultaneous'
			ini_F = randn(size(m_F));
			estTikhonov = ChannelGainMapEstimator('mu_f',1e-4,'ch_reg_f_type','tikhonov','h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'m_tikhonov',m_tikhonov);
			estL1PCO = ChannelGainMapEstimator('mu_f',1e-4,'ch_reg_f_type','l1_PCO','h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType);
			estTotalvariation = ChannelGainMapEstimator('mu_f',1e-5,'ch_reg_f_type','totalvariation','h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','rho',1e-4,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType);
					
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements] = dataGenerator.realization();
			
			% B) estimation
			[m_F_estTikhonov] = estTikhonov.estimate(m_sensorPos,m_sensorInd,v_measurements,[]);
			[m_F_estL1PCO] = estL1PCO.estimate(m_sensorPos,m_sensorInd,v_measurements,[]);
			[m_F_estTotalvariation] = estTotalvariation.estimate(m_sensorPos,m_sensorInd,v_measurements,[]);
			
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_estTikhonov),'tit','Estimated (Tikhonov)');		
			F3 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_estL1PCO),'tit','Estimated (L1/PCO)');
			F4 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_estTotalvariation),'tit','Estimated (TV)');

			F = F_figure('multiplot_array',[F1 F2; F3 F4],'pos',[100.0000  402.0000  [592.3077  592.3077]]);		
			
        end
        
        
		% This is a toy simulation for a non-blind shadow loss field
		% estimation. Independent calibration for Tx/Rx gains and path loss exponent
		% (estimator only knows sensor locations, (shadowing/non-shadowing) channel gain measurements.
		% Non-shadowing measurements are used to estimate the Tx/Rx gains and path loss exponent 
		% vias a least-squares estimator.
		% To run this code, simply set s_pathLossExponent as []
		% for the argument of the estimator "ChannelGainMapEstimator".
		
		function F = compute_fig_2003(obj,niter)
						
			% Create data genator
			m_F = csvread('Map_15_15.csv');
			m_F = m_F/max(max(m_F));
			s_numGrid = size(m_F,1) * size(m_F,2);
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			h_w = @(phi1,phi2) (1/sqrt(phi1))*(phi2<phi1+lambda_W/2); % normalized ellipse model
			s_measurementNum = 300;
            s_pathLossExponent = 2;
            s_sensorNum = 120;
            s_avgSensorGain = 20;
            v_gains = s_avgSensorGain + rand(s_sensorNum,1) .* (2 * (rand(s_sensorNum,1) < 0.5) - 1);
			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
			
			% Create estimator
			ch_reg_f_type = 'l1_PCO'; %'l1_PCO'; %'totalvariation';			
			ch_calibrationType = 'simultaneous'; % 'none','previous','simultaneous'
            m_tikhonov = eye(s_numGrid);
			mu_f = 1e-4; % 1e-4 for l1, 1e-5 for total variation
			ini_F = randn(size(m_F));
            rho =  2.5; % 2.5 for ISTA, 1e-4 for total variation
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'m_tikhonov',m_tikhonov);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
			
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');			
			F = F_figure('multiplot_array',[F1 F2]);
					
        end
		
		
		
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
