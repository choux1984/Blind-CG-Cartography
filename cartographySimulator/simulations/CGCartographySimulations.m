%
%  FIGURES FOR THE PAPER ON BLIND CHANNEL GAIN CARTOGRAPHY
%
%  In this file, we essentially set parameters and invoke a simulator
%

classdef CGCartographySimulations < simFunctionSet
	
	properties
		v_rangPhi1 = [4;8];
		v_rangPhi2 = [4;8];
		v_intervGrid = [0.2;0.01];
		v_Phi1 = 4:0.2:6.3;
		v_Phi2Axis = 4:0.01:6.5;
		
		v_30rangPhi1 = [17;22];
		v_30rangPhi2 = [17;22];
		v_30intervGrid = [0.2;0.01];
		v_30Phi1 = 17:0.2:19.3;
		v_30Phi2Axis = 17:0.01:19.5;
		
		
		v_rangPhi1real = [15;25];
		v_rangPhi2real = [15;27];
		v_intervGridreal = [0.5;0.01];
% 		v_Phi1real = 17.5:0.5:23;
% 		v_Phi2Axisreal = 17.5:0.01:25;	
		
		v_Phi1real = 22:0.5:27.5;
		v_Phi2Axisreal = 22:0.01:29.5;
		
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
		% estimation. The dataset is generated with the normalized ellipse
		% model. Ground truth sensors gains and pathloss exponent are used. 
		function F = compute_fig_2001(obj,niter)
						
			% Create data genator
			m_F = csvread('Map_15_15.csv');
			m_F = m_F/max(max(m_F));
			[s_sizeY,s_sizeX] = size(m_F);
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			s_measurementNum = 200;
            s_pathLossExponent = 2;
            v_gains = [];
			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
			
			% Create estimator
			ch_reg_f_type = 'tikhonov'; %'l1_PCO'; %'totalvariation';
			m_tikhonov = eye(s_sizeX *s_sizeY);
			ch_estimationType = 'non-blind';
			mu_f = 1e-4;
			ini_F = randn(size(m_F));
            rho =  1e-2; % for ISTA, rho is roughtly 2.5
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_estimationType',ch_estimationType,'lambda_W',lambda_W,'m_tikhonov',m_tikhonov);
						
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements);
			
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
			[s_sizeY,s_sizeX] = size(m_F);
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			s_measurementNum = 200;
            s_pathLossExponent = 2;
            v_gains = [];
			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
			
			% Create estimator
			ch_calibrationType = 'none';
			ch_estimationType = 'non-blind';
			ini_F = randn(size(m_F));
			m_tikhonov = eye(s_sizeX *s_sizeY);
			estTikhonov = ChannelGainMapEstimator('mu_f',1e-4,'ch_reg_f_type','tikhonov','h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'ch_estimationType',ch_estimationType,'lambda_W',lambda_W,'m_tikhonov',m_tikhonov);
			estL1PCO = ChannelGainMapEstimator('mu_f',1e-4,'ch_reg_f_type','l1_PCO','h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'ch_estimationType',ch_estimationType,'lambda_W',lambda_W);
			estTotalvariation = ChannelGainMapEstimator('mu_f',1e-5,'ch_reg_f_type','totalvariation','h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','rho',1e-4,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'ch_estimationType',ch_estimationType,'lambda_W',lambda_W);
					
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements] = dataGenerator.realization();
			
			% B) estimation
			[m_F_estTikhonov] = estTikhonov.estimate(m_sensorPos,m_sensorInd,v_measurements);
			[m_F_estL1PCO] = estL1PCO.estimate(m_sensorPos,m_sensorInd,v_measurements);
			[m_F_estTotalvariation] = estTotalvariation.estimate(m_sensorPos,m_sensorInd,v_measurements);
			
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
			[s_sizeY,s_sizeX] = size(m_F);
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
			m_tikhonov = eye(s_sizeX *s_sizeY);
			ch_calibrationType = 'none'; % 'none','simultaneous'
			ch_estimationType = 'non-blind';
			mu_f = 1e-4; % 1e-4 for l1, 1e-5 for total variation
			ini_F = randn(size(m_F));
            rho =  2.5; % 2.5 for ISTA, 1e-4 for total variation
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'ch_estimationType',ch_estimationType,'lambda_W',lambda_W,'m_tikhonov',m_tikhonov);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			est.m_Omega = est.sensorMapOp(m_sensorInd,s_sensorNum);
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(m_sensorPos,m_sensorInd,v_measurementsNoShadowing,est.m_Omega) ;
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements);
			
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');			
			F = F_figure('multiplot_array',[F1 F2]);
					
		end	
		
		% This is a toy simulation of a non-blind shadow loss field
		% estimation with joint sensor gain and pathloss exponent calibration.
		
		function F = compute_fig_2004(obj,niter)
			
			% Create data genator
			m_F = csvread('Map_15_15.csv');
			m_F = m_F/max(max(m_F));
			[s_sizeY,s_sizeX] = size(m_F);
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			s_measurementNum = 400;
			s_pathLossExponent = 2;
			s_sensorNum = 120;
			s_avgSensorGain = 20;
			v_gains = s_avgSensorGain + rand(s_sensorNum,1) .* (2 * (rand(s_sensorNum,1) < 0.5) - 1);
			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
			
			% Create estimator
			ch_reg_f_type = 'l1_PCO'; %'l1_PCO'; %'totalvariation';
			m_tikhonov = eye(s_sizeX *s_sizeY);
			ch_calibrationType = 'simultaneous'; % 'none','simultaneous'
			ch_estimationType = 'non-blind';
			mu_f = 1e-4; % 1e-4 for l1, 1e-5 for total variation
			ini_F = randn(size(m_F));
			rho =  2.5; % 2.5 for ISTA, 1e-4 for total variation
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'ch_estimationType',ch_estimationType,'lambda_W',lambda_W,'m_tikhonov',m_tikhonov);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			est.m_Omega = est.sensorMapOp(m_sensorInd,s_sensorNum);
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(m_sensorPos,m_sensorInd,v_measurementsNoShadowing,est.m_Omega) ;
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements);
			
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');
			F = F_figure('multiplot_array',[F1 F2]);
			
		end
		
		% In this simulation, datasets are generated by using an inverse
		% area ellipse model. 
		function F = compute_fig_2005(obj,niter)
						
			% Create data genator
			m_F = csvread('Map_15_15.csv');
			m_F = m_F/max(max(m_F));
			[s_sizeY,s_sizeX] = size(m_F);
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
			m_tikhonov = eye(s_sizeX *s_sizeY);
			ch_calibrationType = 'none'; % 'none','simultaneous'
			ch_estimationType = 'non-blind';
			mu_f = 1e-4; % 1e-4 for l1, 1e-5 for total variation
			ini_F = randn(size(m_F));
            rho =  2.5; % 2.5 for ISTA, 1e-4 for total variation
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType','non-blind','rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'ch_estimationType',ch_estimationType,'lambda_W',lambda_W,'m_tikhonov',m_tikhonov);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			est.m_Omega = est.sensorMapOp(m_sensorInd,s_sensorNum);
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(m_sensorPos,m_sensorInd,v_measurementsNoShadowing,est.m_Omega) ;
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements);
			
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');			
			F = F_figure('multiplot_array',[F1 F2]);
					
		end		

		% This is a test for spatial loss field reconstruction in higher resolution.
		function F = compute_fig_2006(obj,niter)
						
			% Create data genator
			m_F = csvread('Map_15_15.csv');
			m_F = m_F/max(max(m_F));
			[s_sizeY,s_sizeX] = size(m_F);
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
			s_xAxiSize = (s_sizeX-1)*s_resolution+1;
			s_yAxiSize = (s_sizeY-1)*s_resolution+1;
			ini_F = randn(s_yAxiSize,s_xAxiSize);
			m_tikhonov = eye(s_xAxiSize *s_yAxiSize);
            rho =  1e-3; % for ISTA, rho is roughtly 2.5
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_estimationType',ch_estimationType,'s_resolution',s_resolution,'lambda_W',lambda_W,'m_tikhonov',m_tikhonov);
						
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements);
			
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');			
			F = F_figure('multiplot_array',[F1 F2]);		
			
        end

		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% %%  3. Simple blind simulations with synthetic data
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
		% This is a toy simulation for a blind shadow loss field
		% data: Inverse area ellipse model; estimation only considers grid
		% points within an ellipse.
		% Data calibration required before estimation steps for the SLF and weight function. 
		% Tikhonov reg. is adopted for the estimation of the SLF.
		function F = compute_fig_3001(obj,niter)
			% Create data generator
			m_F = csvread('Map_10_10.csv');
			m_F = m_F/max(max(m_F));
			[s_sizeY,s_sizeX] = size(m_F);
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			%%% Inverse area ellipse model
			s_delta = 3*1e-2;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<=phi1+lambda_W/2);
			%%%
			s_measurementNum = 600;
            s_pathLossExponent = 2;
            s_sensorNum = 400;
            s_avgSensorGain = 20;
            v_gains = s_avgSensorGain + rand(s_sensorNum,1) .* (2 * (rand(s_sensorNum,1) < 0.5) - 1);
			s_noiseVar = 0.0001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
							
			% Create estimator
			s_clusterNum = 800;
			ch_reg_f_type = 'tikhonov'; %'tikhonov','l1_PCO'; %'totalvariation';
			m_tikhonov = eye(s_sizeX *s_sizeY);
			ch_calibrationType = 'none'; % 'none','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_SemiAxisLength4Sample = lambda_W;
			s_kernelStd = 0.08;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(s_kernelStd));
			mu_f = 1e-3; 
			mu_w = 2 * 1e-3;
			ini_F = rand(size(m_F));
            rho =  1e-4; 
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'s_SemiAxisLength4Sample',s_SemiAxisLength4Sample,'m_tikhonov',m_tikhonov);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			% m_Omega needs to be loaded or calculated. Furthermore, 
			est.m_Omega = est.sensorMapOp(m_sensorInd,s_sensorNum);
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(m_sensorPos,m_sensorInd,v_measurementsNoShadowing,est.m_Omega) ;
			[m_F_est,h_w_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements);
			
			% DISPLAY
			% A) spatial loss fields
			s_mu_out = 1;
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocessfunction(m_F,m_F,s_mu_out),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocessfunction(m_F_est,m_F,s_mu_out),'tit','Estimated');	
			F(1) = F_figure('multiplot_array',[F1 F2]);
			
			
			F3 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F4 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');		
			F(2) = F_figure('multiplot_array',[F3 F4]);
			
			% B) weight functions
			s_lengthPhi1 = length(obj.v_Phi1);
			s_lengthPhi2Axis = length(obj.v_Phi2Axis);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1(s_phi1Ind));
			end		
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1,obj.v_rangPhi2,obj.v_intervGrid);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(3) = F_figure('X',obj.v_Phi2Axis,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			
		end
		
		% This is a toy simulation for a blind shadow loss field
		% data: Inverse area ellipse model; estimation only considers grid
		% points within an ellipse. Grid points outside of the ellipse are
		% approximated to those on the boundary of the ellipse.
		% Data calibration is jointly done with the estimation of the field
		% and weight function.
		% Tikhonov reg. is adopted for the estimation of the SLF.
		function F = compute_fig_3002(obj,niter)
			% Create data generator
			m_F = csvread('Map_10_10.csv');
			m_F = m_F/max(max(m_F));
			[s_sizeY,s_sizeX] = size(m_F);
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			%%% Inverse area ellipse model
			s_delta = 3*1e-2;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<=phi1+lambda_W/2);
			%%%
			s_measurementNum = 1000;
            s_pathLossExponent = 2;
            s_sensorNum = 400;
            s_avgSensorGain = 20;
            v_gains = s_avgSensorGain + rand(s_sensorNum,1) .* (2 * (rand(s_sensorNum,1) < 0.5) - 1);
			s_noiseVar = 0.0001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
							
			% Create estimator
			s_clusterNum = 800;
			ch_reg_f_type = 'tikhonov'; %'tikhonov','l1_PCO'; %'totalvariation';
			m_tikhonov = eye(s_sizeX *s_sizeY);
			ch_calibrationType = 'simultaneous'; % 'none','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_SemiAxisLength4Sample = lambda_W;
			s_kernelStd = 0.08;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(s_kernelStd));
			mu_f = 1e-3; 
			mu_w = 2 * 1e-3;
			ini_F = rand(size(m_F));
            rho =  1e-4; 
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'s_SemiAxisLength4Sample',s_SemiAxisLength4Sample,'m_tikhonov',m_tikhonov);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			% m_Omega needs to be loaded or calculated. Furthermore, 
			est.m_Omega = est.sensorMapOp(m_sensorInd,s_sensorNum);
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(m_sensorPos,m_sensorInd,v_measurementsNoShadowing,est.m_Omega) ;
			[m_F_est,h_w_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements);
			
			% DISPLAY
			% A) spatial loss fields
			s_mu_out = 1;
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocessfunction(m_F,m_F,s_mu_out),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocessfunction(m_F_est,m_F,s_mu_out),'tit','Estimated');	
			F(1) = F_figure('multiplot_array',[F1 F2]);
			
			
			F3 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F4 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');		
			F(2) = F_figure('multiplot_array',[F3 F4]);
			
			% B) weight functions
			s_lengthPhi1 = length(obj.v_Phi1);
			s_lengthPhi2Axis = length(obj.v_Phi2Axis);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1(s_phi1Ind));
			end		
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1,obj.v_rangPhi2,obj.v_intervGrid);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(3) = F_figure('X',obj.v_Phi2Axis,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			
		end

		% data: Inverse area ellipse model; estimation only considers grid
		% points within an ellipse. 
		% Data calibration required before estimation steps for the SLF and weight function.
		% l-1 reg. is adopted for the estimation of the SLF.
		function F = compute_fig_3003(obj,niter)
			% Create data generator
			m_F = csvread('Map_10_10.csv');
			m_F = m_F/max(max(m_F));
			[s_sizeY,s_sizeX] = size(m_F);
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			%%% Inverse area ellipse model
			s_delta = 3*1e-2;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<=phi1+lambda_W/2);
			%%%
			s_measurementNum = 600;
            s_pathLossExponent = 2;
            s_sensorNum = 400;
            s_avgSensorGain = 20;
            v_gains = s_avgSensorGain + rand(s_sensorNum,1) .* (2 * (rand(s_sensorNum,1) < 0.5) - 1);
			s_noiseVar = 0.0001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
							
			% Create estimator
			s_clusterNum = 800;
			ch_reg_f_type = 'l1_PCO'; %'tikhonov','l1_PCO'; %'totalvariation';
			m_tikhonov = eye(s_sizeX *s_sizeY);
			ch_calibrationType = 'none'; % 'none','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_SemiAxisLength4Sample = lambda_W;
			s_kernelStd = 0.1;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(s_kernelStd));
			mu_f = 1e-4; 
			mu_w = 2 * 1e-4;
			ini_F = rand(size(m_F));
            rho =  1e-4; 
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'s_SemiAxisLength4Sample',s_SemiAxisLength4Sample,'m_tikhonov',m_tikhonov);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			% m_Omega needs to be loaded or calculated. Furthermore, 
			est.m_Omega = est.sensorMapOp(m_sensorInd,s_sensorNum);
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(m_sensorPos,m_sensorInd,v_measurementsNoShadowing,est.m_Omega) ;
			[m_F_est,h_w_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements);
			
			% DISPLAY
			% A) spatial loss fields
			s_mu_out = 1;
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocessfunction(m_F,m_F,s_mu_out),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocessfunction(m_F_est,m_F,s_mu_out),'tit','Estimated');	
			F(1) = F_figure('multiplot_array',[F1 F2]);
			
			
			F3 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F4 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');		
			F(2) = F_figure('multiplot_array',[F3 F4]);
			
			% B) weight functions
			s_lengthPhi1 = length(obj.v_Phi1);
			s_lengthPhi2Axis = length(obj.v_Phi2Axis);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1(s_phi1Ind));
			end		
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1,obj.v_rangPhi2,obj.v_intervGrid);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(3) = F_figure('X',obj.v_Phi2Axis,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			
		end	

		% data: Inverse area ellipse model; estimation only considers grid
		% points within an ellipse. 
		% Data calibration required before estimation steps for the SLF and weight function.
		% Total variation reg. is adopted for the estimation of the SLF.
		function F = compute_fig_3004(obj,niter)
			% Create data generator
			m_F = csvread('Map_10_10.csv');
			m_F = m_F/max(max(m_F));
			[s_sizeY,s_sizeX] = size(m_F);
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			%%% Inverse area ellipse model
			s_delta = 3*1e-2;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<=phi1+lambda_W/2);
			%%%
			s_measurementNum = 600;
            s_pathLossExponent = 2;
            s_sensorNum = 400;
            s_avgSensorGain = 20;
            v_gains = s_avgSensorGain + rand(s_sensorNum,1) .* (2 * (rand(s_sensorNum,1) < 0.5) - 1);
			s_noiseVar = 0.0001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
							
			% Create estimator
			s_clusterNum = 800;
			ch_reg_f_type = 'totalvariation'; %'tikhonov','l1_PCO'; %'totalvariation';
			m_tikhonov = eye(s_sizeX *s_sizeY);
			ch_calibrationType = 'none'; % 'none','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_SemiAxisLength4Sample = lambda_W;
			s_kernelStd = 0.1;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(s_kernelStd));
			mu_f = 1e-4; 
			mu_w = 2 * 1e-4;
			ini_F = rand(size(m_F));
            rho =  1e-4; 
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'s_SemiAxisLength4Sample',s_SemiAxisLength4Sample,'m_tikhonov',m_tikhonov);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			% m_Omega needs to be loaded or calculated. Furthermore, 
			est.m_Omega = est.sensorMapOp(m_sensorInd,s_sensorNum);
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(m_sensorPos,m_sensorInd,v_measurementsNoShadowing,est.m_Omega) ;
			[m_F_est,h_w_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements);
			
			% DISPLAY
			% A) spatial loss fields
			s_mu_out = 1;
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocessfunction(m_F,m_F,s_mu_out),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocessfunction(m_F_est,m_F,s_mu_out),'tit','Estimated');	
			F(1) = F_figure('multiplot_array',[F1 F2]);
			
			
			F3 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F4 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');		
			F(2) = F_figure('multiplot_array',[F3 F4]);
			
			% B) weight functions
			s_lengthPhi1 = length(obj.v_Phi1);
			s_lengthPhi2Axis = length(obj.v_Phi2Axis);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1(s_phi1Ind));
			end		
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1,obj.v_rangPhi2,obj.v_intervGrid);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(3) = F_figure('X',obj.v_Phi2Axis,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			
		end
		
		% Dimension of the SLF: 30-by-30
		% data: Inverse area ellipse model; estimation only considers grid
		% points within an ellipse.
		% Data calibration required before estimation steps for the SLF and weight function. 
		% Tikhonov reg. is adopted for the estimation of the SLF.
		function F = compute_fig_3005(obj,niter)
			% Create data generator
			m_F = csvread('Map_30_30.csv');
			m_F = m_F/max(max(m_F));
			[s_sizeY,s_sizeX] = size(m_F);
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			%%% Inverse area ellipse model
			s_delta = 3*1e-2;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<=phi1+lambda_W/2);
			%%%
			s_measurementNum = 2500;
            s_pathLossExponent = 2;
            s_sensorNum = 800;
            s_avgSensorGain = 20;
            v_gains = s_avgSensorGain + rand(s_sensorNum,1) .* (2 * (rand(s_sensorNum,1) < 0.5) - 1);
			s_noiseVar = 0.0001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
							
			% Create estimator
			s_clusterNum = 2500;
			ch_reg_f_type = 'tikhonov'; %'tikhonov','l1_PCO'; %'totalvariation';
			m_tikhonov = eye(s_sizeX *s_sizeY);
			ch_calibrationType = 'none'; % 'none','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_SemiAxisLength4Sample = lambda_W;
			s_kernelStd = 0.05;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(s_kernelStd));
			mu_f = 1e-4; 
			mu_w = 2 * 1e-5;
			ini_F = rand(size(m_F));
            rho =  1e-4; 
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'s_SemiAxisLength4Sample',s_SemiAxisLength4Sample,'m_tikhonov',m_tikhonov);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			% m_Omega needs to be loaded or calculated. Furthermore, 
			est.m_Omega = est.sensorMapOp(m_sensorInd,s_sensorNum);
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(m_sensorPos,m_sensorInd,v_measurementsNoShadowing,est.m_Omega) ;
			[m_F_est,h_w_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements);
			
			% DISPLAY
			% A) spatial loss fields
			s_mu_out = 1;
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocessfunction(m_F,m_F,s_mu_out),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocessfunction(m_F_est,m_F,s_mu_out),'tit','Estimated');	
			F(1) = F_figure('multiplot_array',[F1 F2]);
			
			
			F3 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F4 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');		
			F(2) = F_figure('multiplot_array',[F3 F4]);
			
			% B) weight functions
			s_lengthPhi1 = length(obj.v_Phi1);
			s_lengthPhi2Axis = length(obj.v_Phi2Axis);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1(s_phi1Ind));
			end		
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1,obj.v_rangPhi2,obj.v_intervGrid);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(3) = F_figure('X',obj.v_Phi2Axis,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			
		end
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% %%  4. Non-blind simulations with REAL data
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		% This is a realdata simulation for a non-blind shadow loss field
		% estimation. Independent calibration (path loss exponent and sensor 
		% gains are independently esimated with the nonshadowing dataset).
		% Weight matrices are generated by using the normalized ellipse
		% model.
		
		function F = compute_fig_4001(obj,niter)
						
			% Create data genator
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights. e.g. 0.3937 for 1st fresnel zone
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			
			ch_reg_f_type = 'tikhonov'; %'l1_PCO'; %'totalvariation';
			m_spatialCov = ChannelGainMapEstimator.spatialCovMat(s_yAxiSize,s_xAxiSize,s_resolution);
			m_tikhonov = inv(m_spatialCov);
			ch_calibrationType = 'none'; % 'none','simultaneous'
			ch_estimationType = 'non-blind';
			mu_f = 3 * 1e-2; 
			ini_F = randn(s_yAxiSize,s_xAxiSize);
            rho =  1e-4; % for ISTA, rho is roughtly 2.5
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_resolution',s_resolution,'lambda_W',lambda_W,'m_tikhonov',m_tikhonov);
			est.m_Omega = csvread('m_Omega.csv');
			
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(t_sensorPos(:,:,1),t_sensorInd(:,:,1),v_measurementsNoShadowing,est.m_Omega);
			m_sensorPos = t_sensorPos(:,:,2);
			m_sensorInd = t_sensorInd(:,:,2);
			est.m_Omega = est.m_Omega([1:580,601:end],:);
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd(:,[1:580,601:end]),v_measurements);

			F(1) = F_figure('Z',m_F_est,'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			F(2) = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);

		end

		% Real data simulation with simultaneous calibration
		
		function F = compute_fig_4002(obj,niter)
						
			% Create data genator
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 2;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.3; % parameter to determine the threshold for nonzero weights. e.g. 0.3937 for 1st fresnel zone
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			
			ch_reg_f_type = 'tikhonov'; %'l1_PCO'; %'totalvariation';
			m_spatialCov = ChannelGainMapEstimator.spatialCovMat(s_yAxiSize,s_xAxiSize,s_resolution);
			m_tikhonov = inv(m_spatialCov);
			ch_calibrationType = 'simultaneous'; % 'none','simultaneous'
			ch_estimationType = 'non-blind';
			mu_f = 8 * 1e-3; % 'tikhonov' : 5 * 1e-2 / 'totalvariation': 1e-1
			ini_F = randn(s_yAxiSize,s_xAxiSize);
            rho =  1e-4; % for ISTA, rho is roughtly 2.5
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_resolution',s_resolution,'lambda_W',lambda_W,'m_tikhonov',m_tikhonov);
						
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			est.m_Omega = csvread('m_Omega.csv');
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(t_sensorPos(:,:,1),t_sensorInd(:,:,1),v_measurementsNoShadowing,est.m_Omega);
			m_sensorPos = t_sensorPos(:,:,2);
			m_sensorInd = t_sensorInd(:,:,2);
			est.m_Omega = est.m_Omega([1:580,601:end],:);
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd(:,[1:580,601:end]),v_measurements);

			F = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);

		end
		
		% Sensor gains and pathloss exponent are independently estimated by
		% using measurements from free space.
		% Weight matrices are generated by using the inverse area ellipse
		% model.
		function F = compute_fig_4003(obj,niter)
						
			% Create data genator
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights. e.g. 0.3937 for 1st fresnel zone
			%%% Inverse area ellipse model
			s_delta = 1e-3;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<phi1+lambda_W/2);
			%%%
			
			ch_reg_f_type = 'tikhonov'; %'l1_PCO'; %'totalvariation';
			m_spatialCov = ChannelGainMapEstimator.spatialCovMat(s_yAxiSize,s_xAxiSize,s_resolution);
			m_tikhonov = inv(m_spatialCov);
			ch_calibrationType = 'none'; % 'none','simultaneous'
			ch_estimationType = 'non-blind';
			mu_f = 5* 1e-3; % 'tikhonov' : / 'totalvariation': 
			ini_F = randn(s_yAxiSize,s_xAxiSize);
            rho =  1e-4; % for ISTA, rho is roughtly 2.5
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_resolution',s_resolution,'lambda_W',lambda_W,'m_tikhonov',m_tikhonov);
						
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			est.m_Omega = csvread('m_Omega.csv');
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(t_sensorPos(:,:,1),t_sensorInd(:,:,1),v_measurementsNoShadowing,est.m_Omega);
			m_sensorPos = t_sensorPos(:,:,2);
			m_sensorInd = t_sensorInd(:,:,2);
			est.m_Omega = est.m_Omega([1:580,601:end],:);
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd(:,[1:580,601:end]),v_measurements);

			F = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);

		end

		% Sensor gains and pathloss exponent are jointly estimated with the
		% spatial loss field.
		% Weight matrices are generated by using the inverse area ellipse
		% model.
		function F = compute_fig_4004(obj,niter)
						
			% Create data genator
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 2;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.3; % parameter to determine the threshold for nonzero weights. e.g. 0.3937 for 1st fresnel zone
			%%% Inverse area ellipse model
			s_delta = 1e-3;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<phi1+lambda_W/2);
			%%%
			
			ch_reg_f_type = 'tikhonov'; %'l1_PCO'; %'totalvariation';
			m_spatialCov = ChannelGainMapEstimator.spatialCovMat(s_yAxiSize,s_xAxiSize,s_resolution);
			m_tikhonov = inv(m_spatialCov);
			ch_calibrationType = 'simultaneous'; % 'none','simultaneous'
			ch_estimationType = 'non-blind';
			mu_f = 5* 1e-3; % 'tikhonov' : / 'totalvariation': 
			ini_F = randn(s_yAxiSize,s_xAxiSize);
            rho =  1e-4; % for ISTA, rho is roughtly 2.5
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_resolution',s_resolution,'lambda_W',lambda_W,'m_tikhonov',m_tikhonov);
						
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			est.m_Omega = csvread('m_Omega.csv');
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(t_sensorPos(:,:,1),t_sensorInd(:,:,1),v_measurementsNoShadowing,est.m_Omega);
			m_sensorPos = t_sensorPos(:,:,2);
			m_sensorInd = t_sensorInd(:,:,2);
			est.m_Omega = est.m_Omega([1:580,601:end],:);
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd(:,[1:580,601:end]),v_measurements);

			F = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);

		end
	
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% %%  5. Blind simulations with REAL data
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		% This is a realdata simulation for a non-blind shadow loss field
		% estimation. Independent calibration (path loss exponent and sensor 
		% gains are independently esimated with the nonshadowing dataset).
		% Weight matrices are generated by using the normalized ellipse
		% model.
		
		function F = compute_fig_5001(obj,niter)
						
			% Create data genator
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights. e.g. 0.3937 for 1st fresnel zone
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			
			s_clusterNum = 2000;
			ch_reg_f_type = 'tikhonov'; %'tikhonov','l1_PCO'; %'totalvariation';
			m_spatialCov = ChannelGainMapEstimator.spatialCovMat(s_yAxiSize,s_xAxiSize,s_resolution);
			m_tikhonov = inv(m_spatialCov);
			ch_calibrationType = 'none'; % 'none','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_SemiAxisLength4Sample = lambda_W;
			s_kernelStd = 0.12;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(s_kernelStd));
			mu_f = 1e-3; 
			mu_w = 1e-2;
			ini_F = rand(s_yAxiSize,s_xAxiSize);
            rho =  1e-4; 
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'s_SemiAxisLength4Sample',s_SemiAxisLength4Sample,'m_tikhonov',m_tikhonov,'s_resolution',s_resolution);
			est.m_Omega = csvread('m_Omega.csv');
			
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(t_sensorPos(:,:,1),t_sensorInd(:,:,1),v_measurementsNoShadowing,est.m_Omega);
			m_sensorPos = t_sensorPos(:,:,2);
			m_sensorInd = t_sensorInd(:,:,2);
			est.m_Omega = est.m_Omega([1:580,601:end],:);
			[m_F_est,h_w_est] = est.estimate(m_sensorPos,m_sensorInd(:,[1:580,601:end]),v_measurements);
			
			F(1) = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			F(2) = F_figure('Z',m_F_est,'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);

			% B) weight functions
			s_lengthPhi1 = length(obj.v_Phi1real);
			s_lengthPhi2Axis = length(obj.v_Phi2Axisreal);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
			end		
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1real,obj.v_rangPhi2real,obj.v_intervGridreal);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(3) = F_figure('X',obj.v_Phi2Axisreal,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			

		end		

		function F = compute_fig_5002(obj,niter)
			
			% Create data genator
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights. e.g. 0.3937 for 1st fresnel zone
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			
			s_clusterNum = 2500;
			ch_reg_f_type = 'tikhonov'; %'tikhonov','l1_PCO'; %'totalvariation';
			m_spatialCov = ChannelGainMapEstimator.spatialCovMat(s_yAxiSize,s_xAxiSize,s_resolution);
			m_tikhonov = inv(m_spatialCov);
			ch_calibrationType = 'none'; % 'none','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_SemiAxisLength4Sample = lambda_W;
			s_kernelStd = 0.04;
			h_kernel = @(input1,input2) exp(-norm(input1-input2).^2./(2 * s_kernelStd.^2));
			mu_f = 1e-3;
			mu_w = 1e-4;
			ini_F = rand(s_yAxiSize,s_xAxiSize);
			rho =  1e-4;
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'s_SemiAxisLength4Sample',s_SemiAxisLength4Sample,'m_tikhonov',m_tikhonov,'s_resolution',s_resolution);
			est.m_Omega = csvread('m_Omega.csv');
			
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(t_sensorPos(:,:,1),t_sensorInd(:,:,1),v_measurementsNoShadowing,est.m_Omega);
			m_sensorPos = t_sensorPos(:,:,2);
			m_sensorInd = t_sensorInd(:,:,2);
			est.m_Omega = est.m_Omega([1:580,601:end],:);
			[m_F_est,h_w_est] = est.estimate(m_sensorPos,m_sensorInd(:,[1:580,601:end]),v_measurements);
			
			F(1) = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			F(2) = F_figure('Z',m_F_est,'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			
			% B) weight functions
			s_lengthPhi1 = length(obj.v_Phi1real);
			s_lengthPhi2Axis = length(obj.v_Phi2Axisreal);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
			end
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1real,obj.v_rangPhi2real,obj.v_intervGridreal);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(3) = F_figure('X',obj.v_Phi2Axisreal,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			
			
		end
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% %%  6. Blind simulations with REAL data (warm start)
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		% This is a function to estimate optimal coefficient to
		% represent the weight function based on the normalized ellipse model via kernel regression.
		
		function F = compute_fig_6001(obj,niter)
						
			% Create data genator
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights. e.g. 0.3937 for 1st fresnel zone
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			
			s_clusterNum = 1500;
			ch_reg_f_type = 'tikhonov'; %'tikhonov','l1_PCO'; %'totalvariation';
			m_spatialCov = ChannelGainMapEstimator.spatialCovMat(s_yAxiSize,s_xAxiSize,s_resolution);
			m_tikhonov = inv(m_spatialCov);
			ch_calibrationType = 'none'; % 'none','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_SemiAxisLength4Sample = lambda_W;
			s_kernelStd = 0.08;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(s_kernelStd));
			mu_f = 1e-3; 
			mu_w = 1e-4;
			ini_F = rand(s_yAxiSize,s_xAxiSize);
            rho =  1e-4; 
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'s_SemiAxisLength4Sample',s_SemiAxisLength4Sample,'m_tikhonov',m_tikhonov,'s_resolution',s_resolution);
			est.m_Omega = csvread('m_Omega.csv');
						
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,~,~] = dataGenerator.realization();
			
			% B) estimation
			m_sensorPos = t_sensorPos(:,:,2);
			m_sensorInd = t_sensorInd(:,:,2);
			[~,~,v_alpha,h_w_est] = est.optfunctionEstimation(m_sensorPos,m_sensorInd(:,[1:580,601:end]));

			% DISPLAY of weight functions
			s_lengthPhi1 = length(obj.v_Phi1real);
			s_lengthPhi2Axis = length(obj.v_Phi2Axisreal);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
			end
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1real,obj.v_rangPhi2real,obj.v_intervGridreal);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F = F_figure('X',obj.v_Phi2Axisreal,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');

		end

		% This is a function to estimate the spatial loss field via
		% non-blind approach with the estimated weight function via kernel
		% regression.
		function F = compute_fig_6002(obj,niter)
						
			% Create data genator
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights. e.g. 0.3937 for 1st fresnel zone
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			
			s_clusterNum = 1500;
			ch_reg_f_type = 'tikhonov'; %'tikhonov','l1_PCO'; %'totalvariation';
			m_spatialCov = ChannelGainMapEstimator.spatialCovMat(s_yAxiSize,s_xAxiSize,s_resolution);
			m_tikhonov = inv(m_spatialCov);
			ch_calibrationType = 'none'; % 'none','simultaneous'
			ch_estimationType = 'non-blind';
			ch_clustType = 'random';
			s_SemiAxisLength4Sample = lambda_W;
			s_kernelStd = 0.08;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(s_kernelStd));
			mu_f = 1e-3; 
			mu_w = 1e-4;
			ini_F = rand(s_yAxiSize,s_xAxiSize);
            rho =  1e-4; 
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'s_SemiAxisLength4Sample',s_SemiAxisLength4Sample,'m_tikhonov',m_tikhonov,'s_resolution',s_resolution);
			est.m_Omega = csvread('m_Omega.csv');
						
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(t_sensorPos(:,:,1),t_sensorInd(:,:,1),v_measurementsNoShadowing,est.m_Omega);
			m_sensorPos = t_sensorPos(:,:,2);
			m_sensorInd = t_sensorInd(:,:,2);
			est.m_Omega = est.m_Omega([1:580,601:end],:);
			[~,~,~,h_w_est] = est.optfunctionEstimation(m_sensorPos,m_sensorInd(:,[1:580,601:end]));
			est.h_w = h_w_est;
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd(:,[1:580,601:end]),v_measurements);

			F(1) = F_figure('Z',m_F_est,'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			F(2) = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			
			
			% DISPLAY of weight functions
			s_lengthPhi1 = length(obj.v_Phi1real);
			s_lengthPhi2Axis = length(obj.v_Phi2Axisreal);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
			end
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1real,obj.v_rangPhi2real,obj.v_intervGridreal);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(3) = F_figure('X',obj.v_Phi2Axisreal,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');

		end
	
		% This is a function to simulate the blind algorithm
		% where v_alpha (kernel coefficients to represent the weight function)
		% is initialized with the optimal coefficients estimated via kernel regression 
		% of the weight function (normalized ellipse model.)
		% Gaussian kernel is adopted.
		function F = compute_fig_6003(obj,niter)
			
			% Create data genator
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights. e.g. 0.3937 for 1st fresnel zone
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			
			s_clusterNum = 1000;
			ch_reg_f_type = 'tikhonov'; %'tikhonov','l1_PCO'; %'totalvariation';
			m_spatialCov = ChannelGainMapEstimator.spatialCovMat(s_yAxiSize,s_xAxiSize,s_resolution);
			m_tikhonov = inv(m_spatialCov);
			ch_calibrationType = 'none'; % 'none','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_SemiAxisLength4Sample = lambda_W;
			s_kernelStd = 0.001;
			h_kernel = @(input1,input2) exp(-norm(input1-input2).^2./(s_kernelStd));
			mu_f = 1e-3;
			mu_w = 3 * 1e-3;
			ini_F = rand(s_yAxiSize,s_xAxiSize);
			rho =  1e-4;
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'s_SemiAxisLength4Sample',s_SemiAxisLength4Sample,'m_tikhonov',m_tikhonov,'s_resolution',s_resolution);
			est.m_Omega = csvread('m_Omega.csv');
			
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(t_sensorPos(:,:,1),t_sensorInd(:,:,1),v_measurementsNoShadowing,est.m_Omega);
			m_sensorPos = t_sensorPos(:,:,2);
			m_sensorInd = t_sensorInd(:,:,2);
			est.m_Omega = est.m_Omega([1:580,601:end],:);
			[~,~,v_alpha,h_w_est] = est.optfunctionEstimation(m_sensorPos,m_sensorInd(:,[1:580,601:end]));
			est.ini_alpha = v_alpha;
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd(:,[1:580,601:end]),v_measurements);
			
			F(1) = F_figure('Z',m_F_est,'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			F(2) = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			
			
			% DISPLAY of weight functions
			s_lengthPhi1 = length(obj.v_Phi1real);
			s_lengthPhi2Axis = length(obj.v_Phi2Axisreal);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
			end
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1real,obj.v_rangPhi2real,obj.v_intervGridreal);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(3) = F_figure('X',obj.v_Phi2Axisreal,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			
		end

		
		% This is a function to estimate optimal coefficient to
		% represent the weight function based on the normalized ellipse model via kernel regression.
		% Gaussian kernel is adopted to interpolate the weight function
		function F = compute_fig_6004(obj,niter)
			
			% Create data genator
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights. e.g. 0.3937 for 1st fresnel zone
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			
			s_clusterNum = 2500;
			ch_reg_f_type = 'tikhonov'; %'tikhonov','l1_PCO'; %'totalvariation';
			m_spatialCov = ChannelGainMapEstimator.spatialCovMat(s_yAxiSize,s_xAxiSize,s_resolution);
			m_tikhonov = inv(m_spatialCov);
			ch_calibrationType = 'none'; % 'none','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_SemiAxisLength4Sample = lambda_W;
			s_kernelStd = 0.04;
			h_kernel = @(input1,input2) exp(-norm(input1-input2).^2./(2 * s_kernelStd.^2));
			mu_f = 1e-3;
			mu_w = 1e-1;
			ini_F = rand(s_yAxiSize,s_xAxiSize);
			rho =  1e-4;
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'s_SemiAxisLength4Sample',s_SemiAxisLength4Sample,'m_tikhonov',m_tikhonov,'s_resolution',s_resolution);
			est.m_Omega = csvread('m_Omega.csv');
			
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,~,~] = dataGenerator.realization();
			
			% B) estimation
			m_sensorPos = t_sensorPos(:,:,2);
			m_sensorInd = t_sensorInd(:,:,2);
			[~,~,v_alpha,h_w_est] = est.optfunctionEstimation(m_sensorPos,m_sensorInd(:,[1:580,601:end]));
			
			% DISPLAY of weight functions
			s_lengthPhi1 = length(obj.v_Phi1real);
			s_lengthPhi2Axis = length(obj.v_Phi2Axisreal);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
			end
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1real,obj.v_rangPhi2real,obj.v_intervGridreal);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F = F_figure('X',obj.v_Phi2Axisreal,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			
		end
		
		% Initial alpha obatained via kernel regression on normalized ellipse model.
		% Laplacian kernel is adopted.
		% The number of centroids = 700 with only 1 iteration
		% Tikhonov regularization.
		function F = compute_fig_6005(obj,niter)
			
			% Create data genator
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights. e.g. 0.3937 for 1st fresnel zone
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			
			s_clusterNum = 700;
			ch_reg_f_type = 'tikhonov'; %'tikhonov','l1_PCO'; %'totalvariation';
			m_spatialCov = ChannelGainMapEstimator.spatialCovMat(s_yAxiSize,s_xAxiSize,s_resolution);
			m_tikhonov = inv(m_spatialCov);
			ch_calibrationType = 'none'; % 'none','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_SemiAxisLength4Sample = lambda_W;
			s_kernelStd = 0.007;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(s_kernelStd));
			mu_f = 8 * 1e-4;
			mu_w = 7 * 1e-3;
			ini_F = rand(s_yAxiSize,s_xAxiSize);
			rho =  1e-4;
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'s_SemiAxisLength4Sample',s_SemiAxisLength4Sample,'m_tikhonov',m_tikhonov,'s_resolution',s_resolution);
			est.m_Omega = csvread('m_Omega.csv');
			
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(t_sensorPos(:,:,1),t_sensorInd(:,:,1),v_measurementsNoShadowing,est.m_Omega);
			m_sensorPos = t_sensorPos(:,:,2);
			m_sensorInd = t_sensorInd(:,:,2);
			est.m_Omega = est.m_Omega([1:580,601:end],:);
			[~,~,v_alpha,h_w_est] = est.optfunctionEstimation(m_sensorPos,m_sensorInd(:,[1:580,601:end]));
			est.ini_alpha = v_alpha;
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd(:,[1:580,601:end]),v_measurements);
			
			F(1) = F_figure('Z',m_F_est,'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			F(2) = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			
			
			% DISPLAY of weight functions
			s_lengthPhi1 = length(obj.v_Phi1real);
			s_lengthPhi2Axis = length(obj.v_Phi2Axisreal);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
			end
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1real,obj.v_rangPhi2real,obj.v_intervGridreal);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(3) = F_figure('X',obj.v_Phi2Axisreal,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			
		end

		% Initial alpha obatained via kernel regression on normalized ellipse model.
		% Laplacian kernel is adopted.
		% The number of centroids = 5000
		% Tikhonov regularization.
		function F = compute_fig_6006(obj,niter)
			
			% Create data genator
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights. e.g. 0.3937 for 1st fresnel zone
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			
			s_clusterNum = 5000;
			ch_reg_f_type = 'tikhonov'; %'tikhonov','l1_PCO'; %'totalvariation';
			m_spatialCov = ChannelGainMapEstimator.spatialCovMat(s_yAxiSize,s_xAxiSize,s_resolution);
			m_tikhonov = inv(m_spatialCov);
			ch_calibrationType = 'none'; % 'none','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_SemiAxisLength4Sample = lambda_W;
			s_kernelStd = 0.007;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(s_kernelStd));
			mu_f = 4 * 1e-4;
			mu_w = 5 * 1e-3;
			ini_F = rand(s_yAxiSize,s_xAxiSize);
			rho =  1e-4;
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'s_SemiAxisLength4Sample',s_SemiAxisLength4Sample,'m_tikhonov',m_tikhonov,'s_resolution',s_resolution);
			est.m_Omega = csvread('m_Omega.csv');
			
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(t_sensorPos(:,:,1),t_sensorInd(:,:,1),v_measurementsNoShadowing,est.m_Omega);
			m_sensorPos = t_sensorPos(:,:,2);
			m_sensorInd = t_sensorInd(:,:,2);
			est.m_Omega = est.m_Omega([1:580,601:end],:);
			[~,~,v_alpha,h_w_est] = est.optfunctionEstimation(m_sensorPos,m_sensorInd(:,[1:580,601:end]));
			est.ini_alpha = v_alpha;
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd(:,[1:580,601:end]),v_measurements);
			
			F(1) = F_figure('Z',m_F_est,'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			F(2) = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			
			
			% DISPLAY of weight functions
			s_lengthPhi1 = length(obj.v_Phi1real);
			s_lengthPhi2Axis = length(obj.v_Phi2Axisreal);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
			end
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1real,obj.v_rangPhi2real,obj.v_intervGridreal);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(3) = F_figure('X',obj.v_Phi2Axisreal,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			
		end
		
		% Initial alpha obatained via kernel regression on normalized ellipse model.
		% Laplacian kernel is adopted.
		% The number of centroids = 5000
		% Total variation regularization. 
		function F = compute_fig_6007(obj,niter)
			
			% Create data genator
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights. e.g. 0.3937 for 1st fresnel zone
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			
			s_clusterNum = 5000;
			ch_reg_f_type = 'totalvariation'; %'tikhonov','l1_PCO'; %'totalvariation';
			m_spatialCov = ChannelGainMapEstimator.spatialCovMat(s_yAxiSize,s_xAxiSize,s_resolution);
			m_tikhonov = inv(m_spatialCov);
			ch_calibrationType = 'none'; % 'none','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_SemiAxisLength4Sample = lambda_W;
			s_kernelStd = 0.005;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(s_kernelStd));
			mu_f = 1e-3;
			mu_w = 1e-3;
			ini_F = rand(s_yAxiSize,s_xAxiSize);
			rho =  1e-5;
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'s_SemiAxisLength4Sample',s_SemiAxisLength4Sample,'m_tikhonov',m_tikhonov,'s_resolution',s_resolution);
			est.m_Omega = csvread('m_Omega.csv');
			
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(t_sensorPos(:,:,1),t_sensorInd(:,:,1),v_measurementsNoShadowing,est.m_Omega);
			m_sensorPos = t_sensorPos(:,:,2);
			m_sensorInd = t_sensorInd(:,:,2);
			est.m_Omega = est.m_Omega([1:580,601:end],:);
			[~,~,v_alpha,h_w_est] = est.optfunctionEstimation(m_sensorPos,m_sensorInd(:,[1:580,601:end]));
			est.ini_alpha = v_alpha;
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd(:,[1:580,601:end]),v_measurements);
			
			F(1) = F_figure('Z',m_F_est,'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			F(2) = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			
			
			% DISPLAY of weight functions
			s_lengthPhi1 = length(obj.v_Phi1real);
			s_lengthPhi2Axis = length(obj.v_Phi2Axisreal);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
			end
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1real,obj.v_rangPhi2real,obj.v_intervGridreal);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(3) = F_figure('X',obj.v_Phi2Axisreal,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			
		end
		
		% Initial alpha obatained via kernel regression on normalized ellipse model.
		% Laplacian kernel is adopted.
		% The number of centroids = 5000
		% Tikhonov regularization.
		% Pathloss exponent and sensor gains are jointly estimated with
		% v_alpha and v_f.
		function F = compute_fig_6008(obj,niter)
			
			% Create data genator
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights. e.g. 0.3937 for 1st fresnel zone
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			
			s_clusterNum = 5000;
			ch_reg_f_type = 'tikhonov'; %'tikhonov','l1_PCO'; %'totalvariation';
			m_spatialCov = ChannelGainMapEstimator.spatialCovMat(s_yAxiSize,s_xAxiSize,s_resolution);
			m_tikhonov = inv(m_spatialCov);
			ch_calibrationType = 'simultaneous'; % 'none','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_SemiAxisLength4Sample = lambda_W;
			s_kernelStd = 0.005;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(s_kernelStd));
			mu_f = 1e-3;
			mu_w = 1e-3;
			ini_F = rand(s_yAxiSize,s_xAxiSize);
			rho =  1e-4;
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'s_SemiAxisLength4Sample',s_SemiAxisLength4Sample,'m_tikhonov',m_tikhonov,'s_resolution',s_resolution);
			est.m_Omega = csvread('m_Omega.csv');
			
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(t_sensorPos(:,:,1),t_sensorInd(:,:,1),v_measurementsNoShadowing,est.m_Omega);
			m_sensorPos = t_sensorPos(:,:,2);
			m_sensorInd = t_sensorInd(:,:,2);
			est.m_Omega = est.m_Omega([1:580,601:end],:);
			[~,~,v_alpha,h_w_est] = est.optfunctionEstimation(m_sensorPos,m_sensorInd(:,[1:580,601:end]));
			est.ini_alpha = v_alpha;
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd(:,[1:580,601:end]),v_measurements);
			
			F(1) = F_figure('Z',m_F_est,'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			F(2) = F_figure('Z',ChannelGainMapEstimator.postprocessRealJoint(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			
			
			% DISPLAY of weight functions
			s_lengthPhi1 = length(obj.v_Phi1real);
			s_lengthPhi2Axis = length(obj.v_Phi2Axisreal);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
			end
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1real,obj.v_rangPhi2real,obj.v_intervGridreal);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(3) = F_figure('X',obj.v_Phi2Axisreal,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			
		end
		
		% Initial alpha obatained via kernel regression on normalized ellipse model.
		% Laplacian kernel is adopted.
		% The number of centroids = 8000
		% Tikhonov regularization.
		function F = compute_fig_6009(obj,niter)
			
			% Create data genator
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights. e.g. 0.3937 for 1st fresnel zone
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			
			s_clusterNum = 8000;
			ch_reg_f_type = 'tikhonov'; %'tikhonov','l1_PCO'; %'totalvariation';
			m_spatialCov = ChannelGainMapEstimator.spatialCovMat(s_yAxiSize,s_xAxiSize,s_resolution);
			m_tikhonov = inv(m_spatialCov);
			ch_calibrationType = 'none'; % 'none','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_SemiAxisLength4Sample = lambda_W;
			s_kernelStd = 0.01;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(s_kernelStd));
			mu_f = 1e-3;
			mu_w = 4 * 1e-3;
			ini_F = rand(s_yAxiSize,s_xAxiSize);
			rho =  1e-4;
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'s_SemiAxisLength4Sample',s_SemiAxisLength4Sample,'m_tikhonov',m_tikhonov,'s_resolution',s_resolution);
			est.m_Omega = csvread('m_Omega.csv');
			
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(t_sensorPos(:,:,1),t_sensorInd(:,:,1),v_measurementsNoShadowing,est.m_Omega);
			m_sensorPos = t_sensorPos(:,:,2);
			m_sensorInd = t_sensorInd(:,:,2);
			est.m_Omega = est.m_Omega([1:580,601:end],:);
			[~,~,v_alpha,h_w_est] = est.optfunctionEstimation(m_sensorPos,m_sensorInd(:,[1:580,601:end]));
			est.ini_alpha = v_alpha;
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd(:,[1:580,601:end]),v_measurements);
			
			F(1) = F_figure('Z',m_F_est,'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			F(2) = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			
			
			% DISPLAY of weight functions
			s_lengthPhi1 = length(obj.v_Phi1real);
			s_lengthPhi2Axis = length(obj.v_Phi2Axisreal);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
			end
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1real,obj.v_rangPhi2real,obj.v_intervGridreal);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(3) = F_figure('X',obj.v_Phi2Axisreal,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			
		end
		
		
		% This is a function to simulate the blind algorithm
		% where v_alpha (kernel coefficients to represent the weight function)
		% is initialized with the optimal coefficients estimated via kernel regression
		% of the weight function (inverse area ellipse model.)
		% Gaussian kernel is adopted.
		function F = compute_fig_6010(obj,niter)
			
			% Create data genator
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights. e.g. 0.3937 for 1st fresnel zone
			%%% Inverse area ellipse model
			s_delta = 3 * 1e-2;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<phi1+lambda_W/2);
			%%%
						
			s_clusterNum = 5000;
			ch_reg_f_type = 'tikhonov'; %'tikhonov','l1_PCO'; %'totalvariation';
			m_spatialCov = ChannelGainMapEstimator.spatialCovMat(s_yAxiSize,s_xAxiSize,s_resolution);
			m_tikhonov = inv(m_spatialCov);
			ch_calibrationType = 'none'; % 'none','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_SemiAxisLength4Sample = lambda_W;
			s_kernelStd = 0.001;
			h_kernel = @(input1,input2) exp(-norm(input1-input2).^2./(s_kernelStd));
			mu_f = 1e-3;
			mu_w = 1e-3;
			ini_F = rand(s_yAxiSize,s_xAxiSize);
			rho =  1e-4;
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'s_SemiAxisLength4Sample',s_SemiAxisLength4Sample,'m_tikhonov',m_tikhonov,'s_resolution',s_resolution);
			est.m_Omega = csvread('m_Omega.csv');
			
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(t_sensorPos(:,:,1),t_sensorInd(:,:,1),v_measurementsNoShadowing,est.m_Omega);
			m_sensorPos = t_sensorPos(:,:,2);
			m_sensorInd = t_sensorInd(:,:,2);
			est.m_Omega = est.m_Omega([1:580,601:end],:);
			[~,~,v_alpha,h_w_est] = est.optfunctionEstimation(m_sensorPos,m_sensorInd(:,[1:580,601:end]));
			est.ini_alpha = v_alpha;
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd(:,[1:580,601:end]),v_measurements);
			
			F(1) = F_figure('Z',m_F_est,'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			F(2) = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			
			
			% DISPLAY of weight functions
			s_lengthPhi1 = length(obj.v_Phi1real);
			s_lengthPhi2Axis = length(obj.v_Phi2Axisreal);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
			end
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1real,obj.v_rangPhi2real,obj.v_intervGridreal);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(3) = F_figure('X',obj.v_Phi2Axisreal,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			
		end
		
		% Initial alpha obatained via kernel regression on normalized ellipse model.
		% Laplacian kernel is adopted.
		% The number of centroids = 2000
		% Tikhonov regularization.
		% s_resolution = 2(42-by-42 m_est_F)
		function F = compute_fig_6011(obj,niter)
			
			% Create data genator
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 2;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights. e.g. 0.3937 for 1st fresnel zone
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			
			s_clusterNum = 3000;
			ch_reg_f_type = 'tikhonov'; %'tikhonov','l1_PCO'; %'totalvariation';
			m_spatialCov = ChannelGainMapEstimator.spatialCovMat(s_yAxiSize,s_xAxiSize,s_resolution);
			m_tikhonov = inv(m_spatialCov);
			ch_calibrationType = 'none'; % 'none','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_SemiAxisLength4Sample = lambda_W;
			s_kernelStd = 0.01;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(s_kernelStd));
			mu_f = 9 * 1e-4;
			mu_w = 5 * 1e-3;
			ini_F = rand(s_yAxiSize,s_xAxiSize);
			rho =  1e-4;
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'s_SemiAxisLength4Sample',s_SemiAxisLength4Sample,'m_tikhonov',m_tikhonov,'s_resolution',s_resolution);
			est.m_Omega = csvread('m_Omega.csv');
			
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[est.v_gains, est.s_pathLossExponent] = est.estimateSensorGainAndPathLoss(t_sensorPos(:,:,1),t_sensorInd(:,:,1),v_measurementsNoShadowing,est.m_Omega);
			m_sensorPos = t_sensorPos(:,:,2);
			m_sensorInd = t_sensorInd(:,:,2);
			est.m_Omega = est.m_Omega([1:580,601:end],:);
			[~,~,v_alpha,h_w_est] = est.optfunctionEstimation(m_sensorPos,m_sensorInd(:,[1:580,601:end]));
			est.ini_alpha = v_alpha;
			[m_F_est] = est.estimate(m_sensorPos,m_sensorInd(:,[1:580,601:end]),v_measurements);
			
			F(1) = F_figure('Z',m_F_est,'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			F(2) = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			
			
			% DISPLAY of weight functions
			s_lengthPhi1 = length(obj.v_Phi1real);
			s_lengthPhi2Axis = length(obj.v_Phi2Axisreal);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
			end
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1real,obj.v_rangPhi2real,obj.v_intervGridreal);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(3) = F_figure('X',obj.v_Phi2Axisreal,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			
		end
	end
	
	
	
	methods  % to reconsider
		
		function [v_sumOfGains,v_pathLoss] = parameterEstimatorReal(obj,t_sensorPos,t_sensorInd,v_measurementsNoShadowing)
			% Estimate(or assign) the sumOfGain and pathloss for every pair
			% of sensors using structure and non-structured datasets.
			
			% 1. Parameter estimation through non-structured measurements
			m_sensorPosNoShadow=t_sensorPos(:,:,1);
			m_sensorIndNoShadow=t_sensorInd(:,:,1);
			
			s_measurementNum = size(m_sensorIndNoShadow,2);
			v_sensorDistancesNoShadow = zeros(s_measurementNum,1);
			
			for s_measurementInd = 1: s_measurementNum
				v_txPos = m_sensorPosNoShadow(:,m_sensorIndNoShadow(1,s_measurementInd));
				v_rxPos = m_sensorPosNoShadow(:,m_sensorIndNoShadow(2,s_measurementInd));
				v_sensorDistancesNoShadow(s_measurementInd) = 10 * log10(norm(v_txPos-v_rxPos));
			end
			
			m_Omega = csvread('m_Omega_v2.csv');
			s_sensorNum = size(m_Omega,2);
			
			m_regMat = [m_Omega,-1 * v_sensorDistancesNoShadow];
			v_parameters = (m_regMat'*m_regMat + 1e-12* eye(s_sensorNum + 1))\(m_regMat'*v_measurementsNoShadowing);
			v_gains_est = v_parameters(1:s_sensorNum,1);
			s_pathLossExponent_est = v_parameters(s_sensorNum+1,1);
			
			
			
			% 2. determine v_sumOfGains and v_pathLoss for the structured
			% dataset.
			m_sensorPos=t_sensorPos(:,:,2);
			m_sensorInd=t_sensorInd(:,:,2);
			
			v_sensorDistances = zeros(s_measurementNum,1);
			for s_measurementInd = 1: s_measurementNum
				v_txPos = m_sensorPos(:,m_sensorInd(1,s_measurementInd));
				v_rxPos = m_sensorPos(:,m_sensorInd(2,s_measurementInd));
				v_sensorDistances(s_measurementInd) = 10 * log10(norm(v_txPos-v_rxPos));
			end
			
			% 			m_Omega = csvread('m_Omega_v2.csv'); % due to measurements from sensors not conforming to the sensor selction rule
			v_sumOfGains = m_Omega * v_gains_est;
			
			v_pathLoss = s_pathLossExponent_est * v_sensorDistances;
			% Measurements from 581 - 600 are missing in structured test.
			v_sumOfGains = v_sumOfGains([1:580,601:end],1);
			v_pathLoss = v_pathLoss([1:580,601:end],1);
		end
		
	end
	
	
end
