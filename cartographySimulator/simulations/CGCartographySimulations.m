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
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			s_measurementNum = 200;
            s_pathLossExponent = 2;
            v_gains = [];
			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
			
			% Create estimator
			ch_reg_f_type = 'tikhonov'; %'l1_PCO'; %'totalvariation';
			ch_estimationType = 'non-blind';
			mu_f = 1e-4;
			ini_F = randn(size(m_F));
            rho =  1e-2; % for ISTA, rho is roughtly 2.5
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_estimationType',ch_estimationType);
						
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
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			s_measurementNum = 200;
            s_pathLossExponent = 2;
            v_gains = [];
			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
			
			% Create estimator
			ch_calibrationType = 'none'; % 'none','previous','simultaneous'
			ch_estimationType = 'non-blind';
			ini_F = randn(size(m_F));
			estTikhonov = ChannelGainMapEstimator('mu_f',1e-4,'ch_reg_f_type','tikhonov','h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'ch_estimationType',ch_estimationType);
			estL1PCO = ChannelGainMapEstimator('mu_f',1e-4,'ch_reg_f_type','l1_PCO','h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'ch_estimationType',ch_estimationType);
			estTotalvariation = ChannelGainMapEstimator('mu_f',1e-5,'ch_reg_f_type','totalvariation','h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','rho',1e-4,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'ch_estimationType',ch_estimationType);
					
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
		
		function F = compute_fig_2003(obj,niter)
						
			% Create data genator
			m_F = csvread('Map_15_15.csv');
			m_F = m_F/max(max(m_F));
			s_numGrid = size(m_F,1) * size(m_F,2);
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			s_measurementNum = 350;
            s_pathLossExponent = 2;
            s_sensorNum = 120;
            s_avgSensorGain = 20;
            v_gains = s_avgSensorGain + rand(s_sensorNum,1) .* (2 * (rand(s_sensorNum,1) < 0.5) - 1);
			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
			
			% Create estimator
			ch_reg_f_type = 'l1_PCO'; %'l1_PCO'; %'totalvariation';			
			ch_calibrationType = 'simultaneous'; % 'none','previous','simultaneous'
			ch_estimationType = 'non-blind';
			mu_f = 1e-4; % 1e-4 for l1, 1e-5 for total variation
			ini_F = randn(size(m_F));
            rho =  2.5; % 2.5 for ISTA, 1e-4 for total variation
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'ch_estimationType',ch_estimationType);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
			
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');			
			F = F_figure('multiplot_array',[F1 F2]);
					
		end		
		
		% This is a toy example of a non-blind algorithm with an inverse
		% area ellipse model. 
		function F = compute_fig_2004(obj,niter)
						
			% Create data genator
			m_F = csvread('Map_15_15.csv');
			m_F = m_F/max(max(m_F));
			s_numGrid = size(m_F,1) * size(m_F,2);
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			%%% Inverse area ellipse model
			s_delta = 1e-1;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<phi1+lambda_W/2);
			%%%
			s_measurementNum = 300;
            s_pathLossExponent = 2;
            s_sensorNum = 120;
            s_avgSensorGain = 20;
            v_gains = s_avgSensorGain + rand(s_sensorNum,1) .* (2 * (rand(s_sensorNum,1) < 0.5) - 1);
			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
			
			% Create estimator
			ch_reg_f_type = 'l1_PCO'; %'l1_PCO'; %'totalvariation';			
			ch_calibrationType = 'previous'; % 'none','previous','simultaneous'
			ch_estimationType = 'non-blind';
			mu_f = 1e-4; % 1e-4 for l1, 1e-5 for total variation
			ini_F = randn(size(m_F));
            rho =  2.5; % 2.5 for ISTA, 1e-4 for total variation
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'ch_estimationType',ch_estimationType);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
			
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');			
			F = F_figure('multiplot_array',[F1 F2]);
					
		end		

		% Comparison of different regularizers with the inverse area
		% ellipse model. 
		function F = compute_fig_2005(obj,niter)
						
			% Create data genator
			m_F = csvread('Map_15_15.csv');
			m_F = m_F/max(max(m_F));
			s_numGrid = size(m_F,1) * size(m_F,2);
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			%%% Inverse area ellipse model
			s_delta = 1e-1;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<phi1+lambda_W/2);
			%%%
			s_measurementNum = 500;
            s_pathLossExponent = 2;
            s_sensorNum = 120;
            s_avgSensorGain = 20;
            v_gains = s_avgSensorGain + rand(s_sensorNum,1) .* (2 * (rand(s_sensorNum,1) < 0.5) - 1);
			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
			
			% Create estimator
			ch_calibrationType = 'simultaneous'; % 'none','previous','simultaneous'
			ch_estimationType = 'non-blind';
			ini_F = randn(size(m_F));
			estTikhonov = ChannelGainMapEstimator('mu_f',1e-4,'ch_reg_f_type','tikhonov','h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'ch_estimationType',ch_estimationType);
			estL1PCO = ChannelGainMapEstimator('mu_f',1e-4,'ch_reg_f_type','l1_PCO','h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'ch_estimationType',ch_estimationType);
			estTotalvariation = ChannelGainMapEstimator('mu_f',1e-5,'ch_reg_f_type','totalvariation','h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','rho',1e-4,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'ch_estimationType',ch_estimationType);
					
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
		
		
		% This is a test for spatial loss field reconstruction in higher resolution.
		function F = compute_fig_2006(obj,niter)
						
			% Create data genator
			m_F = csvread('Map_15_15.csv');
			m_F = m_F/max(max(m_F));
			[N_x,N_y] = size(m_F);
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			s_measurementNum = 600;
            s_pathLossExponent = 2;
            v_gains = [];
			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
			
			% Create estimator
			ch_reg_f_type = 'tikhonov'; %'l1_PCO'; %'totalvariation';
			ch_estimationType = 'non-blind';
			mu_f = 1e-2;
			s_resolution = 3;
			s_xAxiSize = (N_y-1)*s_resolution+1;
			s_yAxiSize = (N_x-1)*s_resolution+1;
			ini_F = randn(s_yAxiSize,s_xAxiSize);
            rho =  1e-3; % for ISTA, rho is roughtly 2.5
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_estimationType',ch_estimationType,'s_resolution',s_resolution);
						
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements,[]);
			
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');			
			F = F_figure('multiplot_array',[F1 F2]);		
			
        end

		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% %%  3. Simple blind simulations with synthetic data
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
		% This is a toy simulation for a blind shadow loss field
		% estimation with inverse area ellipse model. No need for data
		% calibraion.
		function F = compute_fig_3001(obj,niter)
			% Create data generator
			m_F = csvread('Map_10_10.csv');
			m_F = m_F/max(max(m_F));
			s_numGrid = size(m_F,1) * size(m_F,2);
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			%%% Inverse area ellipse model
			s_delta = 1e-2;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<phi1+lambda_W/2);
			%%%
			s_measurementNum = 1000;
            s_pathLossExponent = 2;
            s_sensorNum = 500;
            s_avgSensorGain = 20;
            v_gains = s_avgSensorGain + rand(s_sensorNum,1) .* (2 * (rand(s_sensorNum,1) < 0.5) - 1);
			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
							
			% Create estimator
			s_clusterNum = 700;
			ch_reg_f_type = 'totalvariation'; %'tikhonov','l1_PCO'; %'totalvariation';			
			ch_calibrationType = 'none'; % 'none','previous','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_kernelStd = 0.2;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(2 * s_kernelStd^2));
			mu_f = 1e-4; % 1e-4 for l1, 1e-5 for total variation
			mu_w = 1e-4; % 1e-4 for l1, 1e-5 for total variation
			ini_F = rand(size(m_F));
            rho =  1e-3; % 2.5 for ISTA, 1e-4 for total variation
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'ch_estimationType',ch_estimationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est,h_w_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
			
			% Display parameter
			v_rangPhi1 = [2;6];
			v_rangPhi2 = [2;6];
			v_intervGrid = [0.2;0.01];
			v_Phi2Axis = 2:0.01:5;
			s_lengthPhi2Axis = length(v_Phi2Axis);
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,v_rangPhi1,v_rangPhi2,v_intervGrid);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');	
			
			F(1) = F_figure('multiplot_array',[F1 F2]);
			F(2) = F_figure('X',v_Phi2Axis,'Y',m_evaluated_w,'colorp',12);

		end
		
		% This is a toy simulation for a blind shadow loss field
		% estimation with normalized ellipse model. No need for data
		% calibraion.
		function F = compute_fig_3002(obj,niter)
			% Create data generator
			m_F = csvread('Map_10_10.csv');
			m_F = m_F/max(max(m_F));
			s_numGrid = size(m_F,1) * size(m_F,2);
			lambda_W = 0.2; % parameter to determine the threshold for nonzero weights
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			s_measurementNum = 1000;
            s_pathLossExponent = 2;
            s_sensorNum = 500;
            s_avgSensorGain = 20;
            v_gains = s_avgSensorGain + rand(s_sensorNum,1) .* (2 * (rand(s_sensorNum,1) < 0.5) - 1);
			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
							
			% Create estimator
			s_clusterNum = 500;
			ch_reg_f_type = 'tikhonov'; %'tikhonov','l1_PCO'; %'totalvariation';			
			ch_calibrationType = 'none'; % 'none','previous','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_kernelStd = 0.18;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(2 * s_kernelStd^2));
			mu_f = 1e-3; % 1e-4 for l1, 1e-5 for total variation
			mu_w = 1e-3; % 1e-4 for l1, 1e-5 for total variation
			ini_F = rand(size(m_F));
            rho =  1e-3; % 2.5 for ISTA, 1e-4 for total variation
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'ch_estimationType',ch_estimationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est,h_w_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
			
			% Display parameter
			v_rangPhi1 = [2;6];
			v_rangPhi2 = [2;6];
			v_intervGrid = [0.2;0.01];	
			[w_o,w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,v_rangPhi1,v_rangPhi2,v_intervGrid);
			
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');	
			
			F = F_figure('multiplot_array',[F1 F2]);
		end
		
		% This is a toy simulation for a blind shadow loss field
		% estimation with inverse area ellipse model. Indepdent data calibration method is used in this simulation.
		
		function F = compute_fig_3003(obj,niter)
			% Create data generator
			m_F = csvread('Map_10_10.csv');
			m_F = m_F/max(max(m_F));
			s_numGrid = size(m_F,1) * size(m_F,2);
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			%%% Inverse area ellipse model
			s_delta = 1e-2;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<phi1+lambda_W/2);
			%%%
			s_measurementNum = 1000;
            s_pathLossExponent = 2;
            s_sensorNum = 500;
            s_avgSensorGain = 20;
            v_gains = s_avgSensorGain + rand(s_sensorNum,1) .* (2 * (rand(s_sensorNum,1) < 0.5) - 1);
			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
							
			% Create estimator
			s_clusterNum = 700;
			ch_reg_f_type = 'totalvariation'; %'tikhonov','l1_PCO'; %'totalvariation';			
			ch_calibrationType = 'previous'; % 'none','previous','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_kernelStd = 0.2;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(2 * s_kernelStd^2));
			mu_f = 1e-4; % 1e-4 for l1, 1e-5 for total variation
			mu_w = 1e-4; % 1e-4 for l1, 1e-5 for total variation
			ini_F = rand(size(m_F));
            rho =  1e-3; % 2.5 for ISTA, 1e-4 for total variation
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'ch_estimationType',ch_estimationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est,h_w_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
			
			% Display parameter
			v_rangPhi1 = [2;6];
			v_rangPhi2 = [2;6];
			v_intervGrid = [0.2;0.01];	
			[w_o,w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,v_rangPhi1,v_rangPhi2,v_intervGrid);
			
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');	
			
			F = F_figure('multiplot_array',[F1 F2]);
% 			F = F_figure('X',(2:0.01:6),'Y',[w_o;w_hat]);

		end
		
		
		% This is a toy simulation for a blind shadow loss field
		% estimation with inverse area ellipse model. Joint data calibration method is used in this simulation.
		function F = compute_fig_3004(obj,niter)
			% Create data generator
			m_F = csvread('Map_10_10.csv');
			m_F = m_F/max(max(m_F));
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			%%% Inverse area ellipse model
			s_delta = 1e-2;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<phi1+lambda_W/2);
			%%%
			s_measurementNum = 1500;
			s_pathLossExponent = 2;
			s_sensorNum = 500;
			s_avgSensorGain = 20;
			v_gains = s_avgSensorGain + rand(s_sensorNum,1) .* (2 * (rand(s_sensorNum,1) < 0.5) - 1);
			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
			
			% Create estimator
			s_clusterNum = 800;
			ch_reg_f_type = 'totalvariation'; %'tikhonov','l1_PCO'; %'totalvariation';
			ch_calibrationType = 'simultaneous'; % 'none','previous','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_kernelStd = 0.2;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(2 * s_kernelStd^2));
			mu_f = 1e-4; % 1e-4 for l1, 1e-5 for total variation
			mu_w = 1e-4; % 1e-4 for l1, 1e-5 for total variation
			ini_F = rand(size(m_F));
			rho =  1e-3; % 2.5 for ISTA, 1e-4 for total variation
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'ch_estimationType',ch_estimationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est,h_w_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
			
			% Display parameter
			v_rangPhi1 = [2;6];
			v_rangPhi2 = [2;6];
			v_intervGrid = [0.2;0.01];
			[w_o,w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,v_rangPhi1,v_rangPhi2,v_intervGrid);
			
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');
			
			F(1) = F_figure('multiplot_array',[F1 F2]);
			% 			F = F_figure('X',(2:0.01:6),'Y',[w_o;w_hat]);
			
			F(2) = F_figure()
		end
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% %%  4. Simple non-blind simulations with REAL data
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		% This is a realdata simulation for a non-blind shadow loss field
		% estimation. Independent calibration (path loss exponent and sensor 
		% gains are independently esimated with the nonshadowing dataset.
		
		function F = compute_fig_4001(obj,niter)
						
			% Create data genator
			ch_dataType = 'real';
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 3;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.3937; % parameter to determine the threshold for nonzero weights
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			
			ch_reg_f_type = 'tikhonov'; %'l1_PCO'; %'totalvariation';
			ch_calibrationType = 'previous'; % should be 'previous','simultaneous'
			ch_estimationType = 'non-blind';
			mu_f = 5 * 1e-2; % 'tikhonov' : 5 * 1e-2 / 'totalvariation': 1e-1
			ini_F = randn(s_yAxiSize,s_xAxiSize);
            rho =  1e-4; % for ISTA, rho is roughtly 2.5
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','rho',rho,'ch_calibrationType',ch_calibrationType,'ch_estimationType',ch_estimationType,'ch_dataType',ch_dataType,'s_resolution',s_resolution);
						
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est] = est.estimate(t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing);			
			
			F = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);

		end
		
		
		% Real data simulation with simultaneous calibration
		
		function F = compute_fig_4002(obj,niter)
						
			% Create data genator
			ch_dataType = 'real';
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 3;
			N_x = 21;
			N_y = 21;
			s_xAxiSize = (N_y-1) * s_resolution + 1;
			s_yAxiSize = (N_x-1) * s_resolution + 1;
			lambda_W = 0.3937; % parameter to determine the threshold for nonzero weights
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			
			ch_reg_f_type = 'tikhonov'; %'l1_PCO'; %'totalvariation';
			ch_calibrationType = 'simultaneous'; % should be 'previous','simultaneous'
			ch_estimationType = 'non-blind';
			mu_f = 5 * 1e-2; % 'tikhonov' : 5 * 1e-2 / 'totalvariation': 1e-1
			ini_F = randn(s_yAxiSize,s_xAxiSize);
            rho =  1e-4; % for ISTA, rho is roughtly 2.5
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','rho',rho,'ch_calibrationType',ch_calibrationType,'ch_estimationType',ch_estimationType,'ch_dataType',ch_dataType,'s_resolution',s_resolution);
						
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est] = est.estimate(t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing);			
			
			F = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);

        end
		
		
		
				
        
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% %%  5. Simple blind simulations with REAL data
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
		
		
		
	end
	
	
end
