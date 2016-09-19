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
		v_Phi1real = 17.5:0.5:23;
		v_Phi2Axisreal = 17.5:0.01:25;	
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
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			%%% Inverse area ellipse model
			s_delta = 2* 1e-2;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<phi1+lambda_W/2);
			%%%
			s_measurementNum = 500;
            s_pathLossExponent = 2;
            s_sensorNum = 500;
            s_avgSensorGain = 20;
            v_gains = s_avgSensorGain + rand(s_sensorNum,1) .* (2 * (rand(s_sensorNum,1) < 0.5) - 1);
			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
							
			% Create estimator
			s_clusterNum = 100;
			ch_reg_f_type = 'tikhonov'; %'tikhonov','l1_PCO'; %'totalvariation';			
			ch_calibrationType = 'none'; % 'none','previous','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_kernelStd = 0.0512;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(s_kernelStd));
			mu_f = 1e-4; % 1e-4 for Tik, 1e-3 for l1, 1e-4 for total variation
			mu_w = 1e-3; % 1e-3 for Tik, 1e-4 for l1, 1e-3 for total variation
			ini_F = rand(size(m_F));
            rho =  1e-4; % 2.5 for ISTA, 1e-4 for total variation
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est,h_w_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
			
			% DISPLAY
			% A) spatial loss fields
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');	
			
			F(1) = F_figure('multiplot_array',[F1 F2]);
			
			% B) weight functions
			s_lengthPhi1 = length(obj.v_Phi1);
			s_lengthPhi2Axis = length(obj.v_Phi2Axis);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1(s_phi1Ind));
			end		
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1,obj.v_rangPhi2,obj.v_intervGrid);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(2) = F_figure('X',obj.v_Phi2Axis,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');

		end
		
		% This is a toy simulation for a blind shadow loss field
		% estimation with normalized ellipse model. No need for data
		% calibraion.
		function F = compute_fig_3002(obj,niter)
			% Create data generator
			m_F = csvread('Map_10_10.csv');
			m_F = m_F/max(max(m_F));
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
			s_clusterNum = 1000;
			ch_reg_f_type = 'totalvariation'; %'tikhonov','l1_PCO'; %'totalvariation';			
			ch_calibrationType = 'none'; % 'none','previous','simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_kernelStd = 0.18;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(2 * s_kernelStd^2));
			mu_f = 1e-4; % 1e-3 for Tikhonov, 1e-4 for l1, 1e-5 for total variation
			mu_w = 1e-3; % 1e-3 for Tikhonov, 1e-4 for l1, 1e-5 for total variation
			ini_F = rand(size(m_F));
            rho =  1e-4; % 2.5 for ISTA, 1e-4 for total variation
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est,h_w_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
			
			% DISPLAY
			% A) spatial loss fields
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');	
			
			F(1) = F_figure('multiplot_array',[F1 F2]);
			
			% B) weight functions
			s_lengthPhi1 = length(obj.v_Phi1);
			s_lengthPhi2Axis = length(obj.v_Phi2Axis);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1(s_phi1Ind));
			end		
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1,obj.v_rangPhi2,obj.v_intervGrid);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(2) = F_figure('X',obj.v_Phi2Axis,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
		end
		
		% This is a toy simulation for a blind shadow loss field
		% estimation with inverse area ellipse model. Independent data calibration method is used in this simulation.
		
		function F = compute_fig_3003(obj,niter)
			% Create data generator
			m_F = csvread('Map_10_10.csv');
			m_F = m_F/max(max(m_F));
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			%%% Inverse area ellipse model
			s_delta = 1e-2;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<phi1+lambda_W/2);
			%%%
			s_measurementNum = 1000;
            s_pathLossExponent = 2;
            s_sensorNum = 300;
            s_avgSensorGain = 20;
%             v_gains = s_avgSensorGain + rand(s_sensorNum,1) .* (2 * (rand(s_sensorNum,1) < 0.5) - 1);
			v_gains = s_avgSensorGain*ones(s_sensorNum,1) + rand(s_sensorNum,1) * 1e-2;

			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
							
			% Create estimator
			s_clusterNum = 400;
			ch_reg_f_type = 'totalvariation'; %'tikhonov'; 'l1_PCO'; 'totalvariation';			
			ch_calibrationType = 'previous'; % 'none';'previous';'simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_kernelStd = 0.2;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(2 * s_kernelStd^2));
			mu_f = 1e-4; % 1e-4 for l1, 1e-5 for total variation
			mu_w = 1e-4; % 1e-4 for l1, 1e-5 for total variation
			ini_F = rand(size(m_F));
            rho =  1e-3; % 2.5 for ISTA, 1e-4 for total variation
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est,h_w_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
					
			% DISPLAY
			% A) spatial loss fields
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');	
			
			F(1) = F_figure('multiplot_array',[F1 F2]);
			
			% B) weight functions
			s_lengthPhi1 = length(obj.v_Phi1);
			s_lengthPhi2Axis = length(obj.v_Phi2Axis);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1(s_phi1Ind));
			end		
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1,obj.v_rangPhi2,obj.v_intervGrid);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(2) = F_figure('X',obj.v_Phi2Axis,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');

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
			s_measurementNum = 2500;
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
			ch_gainType = 'different'; % 'different'; 'same'
			s_kernelStd = 0.2;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(2 * s_kernelStd^2));
			mu_f = 1e-4; % 1e-4 for l1, 1e-5 for total variation
			mu_w = 1e-4; % 1e-4 for l1, 1e-5 for total variation
			ini_F = rand(size(m_F));
			rho =  1e-3; % 2.5 for ISTA, 1e-4 for total variation
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'ch_gainType',ch_gainType);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est,h_w_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
			
			% DISPLAY
			% A) spatial loss fields
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est),'tit','Estimated');	
			
			F(1) = F_figure('multiplot_array',[F1 F2]);
			
			% B) weight functions
			s_lengthPhi1 = length(obj.v_Phi1);
			s_lengthPhi2Axis = length(obj.v_Phi2Axis);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1(s_phi1Ind));
			end		
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1,obj.v_rangPhi2,obj.v_intervGrid);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(2) = F_figure('X',obj.v_Phi2Axis,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
		end
		
		% This is a toy simulation for a blind shadow loss field
		% estimation with inverse area ellipse model. Independent data calibration method is used in this simulation.
		% In addition, the estimated field is displayed after
		% going through a linear postprocessing function.
		function F = compute_fig_3005(obj,niter)
			% Create data generator
			m_F = csvread('Map_15_15.csv');
			m_F = m_F/max(max(m_F));
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			%%% Inverse area ellipse model
			s_delta = 3 * 1e-2;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<phi1+lambda_W/2);
			%%%
			s_measurementNum = 1500;
            s_pathLossExponent = 2;
            s_sensorNum = 300;
            s_avgSensorGain = 20;
%             v_gains = s_avgSensorGain + rand(s_sensorNum,1) .* (2 * (rand(s_sensorNum,1) < 0.5) - 1);
			v_gains = s_avgSensorGain*ones(s_sensorNum,1) + rand(s_sensorNum,1) * 1e-2;

			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
							
			% Create estimator
			s_clusterNum = 1000;
			ch_reg_f_type = 'totalvariation'; %'tikhonov'; 'l1_PCO'; 'totalvariation';			
			ch_calibrationType = 'previous'; % 'none';'previous';'simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_kernelStd = 0.2;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(2 * s_kernelStd^2));
			mu_f = 1e-3; % 1e-4 for l1, 1e-5 for total variation
			mu_w = 1e-4; % 1e-4 for l1, 1e-5 for total variation
			ini_F = rand(size(m_F));
            rho =  1e-4; % 2.5 for ISTA, 1e-4 for total variation
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est,h_w_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
					
			% DISPLAY
			% A) spatial loss fields
			s_mu_out = 1;
			F1 = F_figure('Z',ChannelGainMapEstimator.postprocessfunction(m_F,m_F,s_mu_out),'tit','Original');
			F2 = F_figure('Z',ChannelGainMapEstimator.postprocessfunction(m_F_est,m_F,s_mu_out),'tit','Estimated');	
			
			F(1) = F_figure('multiplot_array',[F1 F2]);
			
			% B) weight functions
			s_lengthPhi1 = length(obj.v_Phi1);
			s_lengthPhi2Axis = length(obj.v_Phi2Axis);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1(s_phi1Ind));
			end		
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1,obj.v_rangPhi2,obj.v_intervGrid);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(2) = F_figure('X',obj.v_Phi2Axis,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');

		end
		
		% This is a toy simulation for a blind shadow loss field
		% data: Inverse area ellipse model; estimation only considers grid
		% points within an ellipse
		% No need for data calibraion.
		
		function F = compute_fig_3006(obj,niter)
			% Create data generator
			m_F = csvread('Map_10_10.csv');
			m_F = m_F/max(max(m_F));
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			%%% Inverse area ellipse model
			s_delta = 3*1e-2;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<phi1+lambda_W/2);
			%%%

			s_measurementNum = 700;
            s_pathLossExponent = 2;
            s_sensorNum = 400;
            s_avgSensorGain = 20;
            v_gains = s_avgSensorGain + rand(s_sensorNum,1) .* (2 * (rand(s_sensorNum,1) < 0.5) - 1);
			s_noiseVar = 0.0001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
							
			% Create estimator
			s_clusterNum = 250;
			ch_reg_f_type = 'tikhonov'; %'tikhonov','l1_PCO'; %'totalvariation';			
			ch_calibrationType = 'none'; % 'none','previous','simultaneous'
			ch_estimationType = 'blind-only-ellipse';
			ch_clustType = 'random';
			s_SemiAxisLength4Sample = lambda_W;
			s_kernelStd = 0.05;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(s_kernelStd));
			mu_f = 1e-4; % 1e-4 for l1, 1e-5 for total variation
			mu_w = 1e-4; % 1e-4 for l1, 1e-5 for total variation
			ini_F = rand(size(m_F));
            rho =  1e-4; % 2.5 for ISTA, 1e-4 for total variation
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'s_SemiAxisLength4Sample',s_SemiAxisLength4Sample);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est,h_w_est] = est.estimate(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
			
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
		
		% This is a synthetic data simulation for journal works with TV regularizer.
		% Channel gain coefficients are estimated from the non-shadowing
		% channel gain measurements. In addition, the inverse area model
		% was adopted.
		function F = compute_fig_3007(obj,niter)
			% Create data generator
			m_F = csvread('Map_30_30.csv');
			m_F = m_F/max(max(m_F));
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			%%% Inverse area ellipse model
			s_delta = 4 * 1e-2;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<phi1+lambda_W/2);
			%%%
			s_measurementNum = 2500;
            s_pathLossExponent = 2;
            s_sensorNum = 600;
            s_avgSensorGain = 20;
%             v_gains = s_avgSensorGain + rand(s_sensorNum,1) .* (2 * (rand(s_sensorNum,1) < 0.5) - 1);
			v_gains = s_avgSensorGain*ones(s_sensorNum,1) + rand(s_sensorNum,1) * 1e-2;

			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
							
			% Create estimator
			s_clusterNum1500 = 1500;
			s_clusterNum2500 = 2500;
			ch_reg_f_type = 'totalvariation'; %'tikhonov'; 'l1_PCO'; 'totalvariation';			
			ch_calibrationType = 'previous'; % 'none';'previous';'simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_kernelStd = 0.0512;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(s_kernelStd));
			mu_f = 1e-4; % 1e-4 for l1, 1e-5 for total variation
			mu_w = 1e-3; % 1e-4 for l1, 1e-5 for total variation
			ini_F = rand(size(m_F));
            rho =  1e-3; % 2.5 for ISTA, 1e-4 for total variation
			est1500 = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum1500,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel);
			est2500 = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum2500,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est1500,h_w_est1500] = est1500.estimate(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
			[m_F_est2500,h_w_est2500] = est2500.estimate(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
					
			% DISPLAY
			% A) spatial loss fields
			F(1) = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est1500),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);	
			F(2) = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est2500),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);	
						
			% B) weight functions
			s_lengthPhi1 = length(obj.v_30Phi1);
			s_lengthPhi2Axis = length(obj.v_30Phi2Axis);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_30Phi1(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_30Phi1(s_phi1Ind));
			end		
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est1500,obj.v_30rangPhi1,obj.v_30rangPhi2,obj.v_30intervGrid);
			m_evaluated_w1500 = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est2500,obj.v_30rangPhi1,obj.v_30rangPhi2,obj.v_30intervGrid);
			m_evaluated_w2500 = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(3) = F_figure('X',obj.v_30Phi2Axis,'Y',m_evaluated_w1500,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			F(4) = F_figure('X',obj.v_30Phi2Axis,'Y',m_evaluated_w2500,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');

		end
		
		% This is a synthetic data simulation for journal works with Tikhonov regularizer.
		% Channel gain coefficients are estimated from the non-shadowing
		% channel gain measurements. In addition, the inverse area model
		% was adopted.
		function F = compute_fig_3008(obj,niter)
			% Create data generator
			m_F = csvread('Map_30_30.csv');
			m_F = m_F/max(max(m_F));
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			%%% Inverse area ellipse model
			s_delta = 4 * 1e-2;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<phi1+lambda_W/2);
			%%%
			s_measurementNum = 2500;
            s_pathLossExponent = 2;
            s_sensorNum = 600;
            s_avgSensorGain = 20;
%             v_gains = s_avgSensorGain + rand(s_sensorNum,1) .* (2 * (rand(s_sensorNum,1) < 0.5) - 1);
			v_gains = s_avgSensorGain*ones(s_sensorNum,1) + rand(s_sensorNum,1) * 1e-2;

			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
							
			% Create estimator
			s_clusterNum1500 = 1500;
			s_clusterNum2500 = 2500;
			ch_reg_f_type = 'tikhonov'; %'tikhonov'; 'l1_PCO'; 'totalvariation';			
			ch_calibrationType = 'previous'; % 'none';'previous';'simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_kernelStd = 0.0512;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(s_kernelStd));
			mu_f = 1e-2; % 1e-4 for l1, 1e-5 for total variation
			mu_w = 1e-1; % 1e-4 for l1, 1e-5 for total variation
			ini_F = rand(size(m_F));
            rho =  1e-3; % 2.5 for ISTA, 1e-4 for total variation
			est1500 = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum1500,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel);
			est2500 = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum2500,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est1500,h_w_est1500] = est1500.estimate(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
			[m_F_est2500,h_w_est2500] = est2500.estimate(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
					
			% DISPLAY
			% A) spatial loss fields
			F(1) = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est1500),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);	
			F(2) = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est2500),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);	
						
			% B) weight functions
			s_lengthPhi1 = length(obj.v_30Phi1);
			s_lengthPhi2Axis = length(obj.v_30Phi2Axis);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_30Phi1(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_30Phi1(s_phi1Ind));
			end		
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est1500,obj.v_30rangPhi1,obj.v_30rangPhi2,obj.v_30intervGrid);
			m_evaluated_w1500 = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est2500,obj.v_30rangPhi1,obj.v_30rangPhi2,obj.v_30intervGrid);
			m_evaluated_w2500 = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(3) = F_figure('X',obj.v_30Phi2Axis,'Y',m_evaluated_w1500,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			F(4) = F_figure('X',obj.v_30Phi2Axis,'Y',m_evaluated_w2500,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');

		end

		% This is a synthetic data simulation for journal works with l-1 regularizer.
		% Channel gain coefficients are estimated from the non-shadowing
		% channel gain measurements. In addition, the inverse area model
		% was adopted.
		function F = compute_fig_3009(obj,niter)
			% Create data generator
			m_F = csvread('Map_30_30.csv');
			m_F = m_F/max(max(m_F));
			lambda_W = 0.4; % parameter to determine the threshold for nonzero weights
			%%% Inverse area ellipse model
			s_delta = 4 * 1e-2;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<phi1+lambda_W/2);
			%%%
			s_measurementNum = 2500;
            s_pathLossExponent = 2;
            s_sensorNum = 600;
            s_avgSensorGain = 20;
%             v_gains = s_avgSensorGain + rand(s_sensorNum,1) .* (2 * (rand(s_sensorNum,1) < 0.5) - 1);
			v_gains = s_avgSensorGain*ones(s_sensorNum,1) + rand(s_sensorNum,1) * 1e-2;

			s_noiseVar = 0.001;% noise variance
			dataGenerator = SyntheticSensorMeasurementsGenerator('m_F',m_F,'h_w',h_w,'s_measurementNum',s_measurementNum,'s_noiseVar',s_noiseVar,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent);
							
			% Create estimator
			s_clusterNum1500 = 1500;
			s_clusterNum2500 = 2500;
			ch_reg_f_type = 'l1_PCO'; %'tikhonov'; 'l1_PCO'; 'totalvariation';			
			ch_calibrationType = 'previous'; % 'none';'previous';'simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_kernelStd = 0.0512;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(s_kernelStd));
			mu_f = 1e-3; % 1e-4 for l1, 1e-5 for total variation
			mu_w = 1e-3; % 1e-4 for l1, 1e-5 for total variation
			ini_F = rand(size(m_F));
            rho =  1.5; % 2.5 for ISTA, 1e-4 for total variation
			est1500 = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum1500,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel);
			est2500 = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'v_gains',v_gains,'s_pathLossExponent',s_pathLossExponent,'ch_calibrationType',ch_calibrationType,'s_clusterNum',s_clusterNum2500,'lambda_W',lambda_W,'ch_clustType',ch_clustType,'h_kernel',h_kernel);
			
			% SIMULATION
			% A) data generation
			[m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est1500,h_w_est1500] = est1500.estimate(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
			[m_F_est2500,h_w_est2500] = est2500.estimate(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
					
			% DISPLAY
			% A) spatial loss fields
			F(1) = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est1500),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);	
			F(2) = F_figure('Z',ChannelGainMapEstimator.postprocess(m_F_est2500),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);	
						
			% B) weight functions
			s_lengthPhi1 = length(obj.v_30Phi1);
			s_lengthPhi2Axis = length(obj.v_30Phi2Axis);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_30Phi1(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_30Phi1(s_phi1Ind));
			end		
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est1500,obj.v_30rangPhi1,obj.v_30rangPhi2,obj.v_30intervGrid);
			m_evaluated_w1500 = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est2500,obj.v_30rangPhi1,obj.v_30rangPhi2,obj.v_30intervGrid);
			m_evaluated_w2500 = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(3) = F_figure('X',obj.v_30Phi2Axis,'Y',m_evaluated_w1500,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			F(4) = F_figure('X',obj.v_30Phi2Axis,'Y',m_evaluated_w2500,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');

		end
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% %%  4. Non-blind simulations with REAL data
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		% This is a realdata simulation for a non-blind shadow loss field
		% estimation. Independent calibration (path loss exponent and sensor 
		% gains are independently esimated with the nonshadowing dataset.
		
		function F = compute_fig_4001(obj,niter)
						
			% Create data genator
			ch_dataType = 'real';
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.3; % parameter to determine the threshold for nonzero weights. e.g. 0.3937 for 1st fresnel zone
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			
			ch_reg_f_type = 'tikhonov'; %'l1_PCO'; %'totalvariation';
			ch_calibrationType = 'previous'; % should be 'previous','simultaneous'
			ch_estimationType = 'non-blind';
			mu_f = 8 * 1e-3; % 'tikhonov' : 5 * 1e-2 / 'totalvariation': 1e-1
			ini_F = randn(s_yAxiSize,s_xAxiSize);
            rho =  1e-4; % for ISTA, rho is roughtly 2.5
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'ch_dataType',ch_dataType,'s_resolution',s_resolution,'lambda_W',lambda_W);
						
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
			s_resolution = 1.5;
			N_x = 21;
			N_y = 21;
			s_xAxiSize = (N_y-1) * s_resolution + 1;
			s_yAxiSize = (N_x-1) * s_resolution + 1;
			lambda_W = 0.3937; % parameter to determine the threshold for nonzero weights
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			
			ch_reg_f_type = 'tikhonov'; %'l1_PCO'; %'totalvariation';
			ch_calibrationType = 'simultaneous'; % should be 'previous','simultaneous'
			ch_estimationType = 'non-blind';
			mu_f = 8 * 1e-3; % 'tikhonov' : 5 * 1e-2 / 'totalvariation': 1e-1
			ini_F = randn(s_yAxiSize,s_xAxiSize);
            rho =  1e-4; % for ISTA, rho is roughtly 2.5
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'ch_dataType',ch_dataType,'s_resolution',s_resolution,'lambda_W',lambda_W);
						
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est] = est.estimate(t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing);			
			
			F = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);

		end
		
		% Real data simulation with the inverse area ellipse model 
		
		function F = compute_fig_4003(obj,niter)
						
			% Create data genator
			ch_dataType = 'real';
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.3937; % parameter to determine the threshold for nonzero weights
			%%% Inverse area ellipse model
			s_delta = 1e-3;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<phi1+lambda_W/2);
			%%%
			
			ch_reg_f_type = 'tikhonov'; %'l1_PCO'; %'totalvariation';
			ch_calibrationType = 'previous'; % should be 'previous','simultaneous'
			ch_estimationType = 'non-blind';
			mu_f = 8 * 1e-3; % 'tikhonov' : 5 * 1e-2 / 'totalvariation': 1e-1
			ini_F = randn(s_yAxiSize,s_xAxiSize);
            rho =  1e-4; % for ISTA, rho is roughtly 2.5
			est = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'ch_dataType',ch_dataType,'s_resolution',s_resolution);
						
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est] = est.estimate(t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing);			
			
			F = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);

		end
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% %%  5. Blind simulations with REAL data
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		function F = compute_fig_5001(obj,niter)
						
			% Create data genator
			ch_dataType = 'real';
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.3937; % parameter to determine the threshold for nonzero weights
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			
			s_clusterNum = 1500;
			ch_reg_f_type = 'tikhonov'; %'tikhonov'; 'l1_PCO'; 'totalvariation';
			ch_calibrationType = 'previous'; % 'none';'previous';'simultaneous'
			ch_estimationType = 'blind';
			ch_clustType = 'random';
			s_kernelStd = 0.09;
			h_kernel = @(input1,input2) exp(-norm(input1-input2).^2./(2 * s_kernelStd^2));
			mu_f = 8 * 1e-3;
			mu_w = 1e-5; % 1e-4 for l1, 1e-5 for total variation
			rho =  1e-4; % for ISTA, rho is roughtly 2.5
			ini_F = rand(s_yAxiSize,s_xAxiSize);
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'ch_dataType',ch_dataType,'s_resolution',s_resolution,'s_clusterNum',s_clusterNum,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'lambda_W',lambda_W);
						
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est,h_w_est] = est.estimate(t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing);
			
			% DISPLAY
			% A) spatial loss fields

			F(1) = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			
			% B) weight functions
			s_lengthPhi1 = length(obj.v_Phi1real);
			s_lengthPhi2Axis = length(obj.v_Phi2Axisreal);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
			end
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1real,obj.v_rangPhi2real,obj.v_intervGridreal);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(2) = F_figure('X',obj.v_Phi2Axisreal,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');

			
		end
		
		
		% This is a function to test the blind algorithm with data after
		% obataining optimal alpha from the algorithm 6002 with the
		% constant weight function
		function F = compute_fig_5002(obj,niter)
						
			% Create data genator
			ch_dataType = 'real';
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.3937; % parameter to determine the threshold for nonzero weights
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			
			s_clusterNum = 2500;
			ch_reg_f_type = 'tikhonov'; %'tikhonov'; 'l1_PCO'; 'totalvariation';
			ch_calibrationType = 'previous'; % 'none';'previous';'simultaneous'
			ch_estimationType = 'warm';
			ch_clustType = 'random';
			ch_coefficient = 'centroids_estimated_constant.mat';
			s_kernelStd = 0.07;
			h_kernel = @(input1,input2) exp(-norm(input1-input2).^2./(2 * s_kernelStd^2));
			mu_f = 1e-4;
			mu_w = 1e-3; % 1e-4 for l1, 1e-5 for total variation
			rho =  1e-4; % for ISTA, rho is roughtly 2.5
			ini_F = rand(s_yAxiSize,s_xAxiSize);
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'ch_dataType',ch_dataType,'s_resolution',s_resolution,'s_clusterNum',s_clusterNum,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'lambda_W',lambda_W,'ch_coefficient',ch_coefficient);
						
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est,h_w_est] = est.estimate(t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing);
			
			% DISPLAY
			% A) spatial loss fields

			F(1) = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			
			% B) weight functions
			s_lengthPhi1 = length(obj.v_Phi1real);
			s_lengthPhi2Axis = length(obj.v_Phi2Axisreal);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
			end
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1real,obj.v_rangPhi2real,obj.v_intervGridreal);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(2) = F_figure('X',obj.v_Phi2Axisreal,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');

			
		end

		% This is a function to test the blind algorithm with data after
		% obataining optimal alpha from the algorithm 6001 with inverse
		% area model
		function F = compute_fig_5003(obj,niter)
						
			% Create data genator
			ch_dataType = 'real';
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.3937; % parameter to determine the threshold for nonzero weights
			%%% Inverse area ellipse model
			s_delta = 1e-3;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<phi1+lambda_W/2);
			%%%
			
			s_clusterNum = 1500;
			ch_reg_f_type = 'tikhonov'; %'tikhonov'; 'l1_PCO'; 'totalvariation';
			ch_calibrationType = 'previous'; % 'none';'previous';'simultaneous'
			ch_estimationType = 'warm';
			ch_clustType = 'random';
			ch_coefficient = 'centroids_estimated_inverse.mat';
			s_kernelStd = 0.05;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(2 * s_kernelStd^2));
			mu_f = 1e-3;
			mu_w = 1e-5; % 1e-4 for l1, 1e-5 for total variation
			rho =  1e-4; % for ISTA, rho is roughtly 2.5
			ini_F = rand(s_yAxiSize,s_xAxiSize);
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'ch_dataType',ch_dataType,'s_resolution',s_resolution,'s_clusterNum',s_clusterNum,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'lambda_W',lambda_W,'ch_coefficient',ch_coefficient);
						
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est,h_w_est] = est.estimate(t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing);
			
			% DISPLAY
			% A) spatial loss fields

			F(1) = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			
			% B) weight functions
			s_lengthPhi1 = length(obj.v_Phi1real);
			s_lengthPhi2Axis = length(obj.v_Phi2Axisreal);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
			end
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1real,obj.v_rangPhi2real,obj.v_intervGridreal);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(2) = F_figure('X',obj.v_Phi2Axisreal,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');

			
		end		
		
		function F = compute_fig_5004(obj,niter)
						
			% Create data genator
			ch_dataType = 'real';
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.3937; % parameter to determine the threshold for nonzero weights
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			
			s_clusterNum = 50;
			ch_reg_f_type = 'tikhonov'; %'tikhonov'; 'l1_PCO'; 'totalvariation';
			ch_calibrationType = 'previous'; % 'none';'previous';'simultaneous'
			ch_estimationType = 'blind-only-ellipse';
			ch_clustType = 'random';
			s_kernelStd = 0.09;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(2 * s_kernelStd^2));
			mu_f = 8 * 1e-3;
			mu_w = 1e-3; % 1e-4 for l1, 1e-5 for total variation
			rho =  1e-4; % for ISTA, rho is roughtly 2.5
			ini_F = randn(s_yAxiSize,s_xAxiSize);
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'rho',rho,'ch_calibrationType',ch_calibrationType,'ch_dataType',ch_dataType,'s_resolution',s_resolution,'s_clusterNum',s_clusterNum,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'lambda_W',lambda_W);
						
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est,h_w_est] = est.estimate(t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing);
			
			% DISPLAY
			% A) spatial loss fields

			F(1) = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
			
			% B) weight functions
			s_lengthPhi1 = length(obj.v_Phi1real);
			s_lengthPhi2Axis = length(obj.v_Phi2Axisreal);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
			end
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1real,obj.v_rangPhi2real,obj.v_intervGridreal);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(2) = F_figure('X',obj.v_Phi2Axisreal,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');

			
		end
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% %%  6. Weight function estimation under RKHS with REAL data 
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		% This is a function to estimate optimal coefficient to
		% represent the weight function of the inverse area ellipse model.
		
		function F = compute_fig_6001(obj,niter)
						
			% Create data genator
			ch_dataType = 'real';
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.3937; % parameter to determine the threshold for nonzero weights 0.3937
			%%% Inverse area ellipse model
			s_delta = 1e-3;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<phi1+lambda_W/2);
			%%%
			s_clusterNum = 1500;
			ch_clustType = 'random';
			ch_estimationType = 'blind';
			s_kernelStd = 0.11;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(2 * s_kernelStd^2));
			mu_w = 1e-5; % 1e-4 for l1, 1e-5 for total variation
			ini_F = rand(s_yAxiSize,s_xAxiSize);
			est = ChannelGainMapEstimator('mu_w',mu_w,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'ch_dataType',ch_dataType,'s_resolution',s_resolution,'s_clusterNum',s_clusterNum,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'lambda_W',lambda_W);
						
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,~,~] = dataGenerator.realization();
			
			% B) estimation
			[~,~,~,h_w_est] = est.optfunctionEstimation(t_sensorPos,t_sensorInd);

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

		% This is a function to estimate optimal coefficient to
		% represent the weight function of the normalized ellipse model.
		function F = compute_fig_6002(obj,niter)
						
			% Create data genator
			ch_dataType = 'real';
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.3937; % parameter to determine the threshold for nonzero weights
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			s_clusterNum = 2500;
			ch_clustType = 'random';
			ch_estimationType = 'blind';
			s_kernelStd = 0.07;
			h_kernel = @(input1,input2) exp(-norm(input1-input2).^2./(2 * s_kernelStd^2));

			mu_f = 8 * 1e-3;
			mu_w = 1e-3; % 1e-4 for l1, 1e-5 for total variation
			rho =  1e-4; % for ISTA, rho is roughtly 2.5
			ini_F = rand(s_yAxiSize,s_xAxiSize);
			est = ChannelGainMapEstimator('mu_w',mu_w,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'ch_dataType',ch_dataType,'s_resolution',s_resolution,'s_clusterNum',s_clusterNum,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'lambda_W',lambda_W);
						
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,~,~] = dataGenerator.realization();
			
			% B) estimation
			[~,~,~,h_w_est] = est.optfunctionEstimation(t_sensorPos,t_sensorInd);

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
		
		% This is a function to estimate optimal coefficient to
		% represent the weight function of the normalized ellipse model. Then, estimate the SLF using the
		% estimated weight function in non-blind fashion.
		
		function F = compute_fig_6003(obj,niter)
						
			% Create data genator
			ch_dataType = 'real';
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.3937; % parameter to determine the threshold for nonzero weights
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			ch_reg_f_type = 'tikhonov'; %'tikhonov'; 'l1_PCO'; 'totalvariation';
			ch_calibrationType = 'previous'; % 'none';'previous';'simultaneous'
			s_clusterNum = 1500;
			ch_clustType = 'random';
			ch_estimationType = 'blind';
			s_kernelStd = 0.11;
			h_kernel = @(input1,input2) exp(-norm(input1-input2).^2./(2 * s_kernelStd^2));

			mu_f = 8 * 1e-3;
			mu_w = 1e-5; % 1e-4 for l1, 1e-5 for total variation
			rho =  1e-4; % for ISTA, rho is roughtly 2.5
			ini_F = rand(s_yAxiSize,s_xAxiSize);
			est = ChannelGainMapEstimator('mu_w',mu_w,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'ch_dataType',ch_dataType,'s_resolution',s_resolution,'s_clusterNum',s_clusterNum,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'lambda_W',lambda_W);
						
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[~,~,~,h_w_est] = est.optfunctionEstimation(t_sensorPos,t_sensorInd);
			est1 = ChannelGainMapEstimator('mu_f',mu_f,'ch_reg_f_type',ch_reg_f_type,'h_w',h_w_est,'ini_F',ini_F,'ch_estimationType','non-blind','rho',rho,'ch_calibrationType',ch_calibrationType,'ch_dataType',ch_dataType,'s_resolution',s_resolution,'lambda_W',lambda_W);
			[m_F_est] = est1.estimate(t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing);			
% 
% 
			% DISPLAY of weight functions
			s_lengthPhi1 = length(obj.v_Phi1real);
			s_lengthPhi2Axis = length(obj.v_Phi2Axisreal);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
			end
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1real,obj.v_rangPhi2real,obj.v_intervGridreal);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(1) = F_figure('X',obj.v_Phi2Axisreal,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			F(2) = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);
	
		end
		
		% This is a function to see changes in fitting on the inverse area
		% weight function with the estimated weight function as the number
		% of cetroids increases.
		function F = compute_fig_6004(obj,niter)
						
			% Create data genator
			ch_dataType = 'real';
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.3937; % parameter to determine the threshold for nonzero weights
			%%% Inverse area ellipse model
			s_delta = 1e-3;
			h_Omega = @(phi1,phi2)  4 ./ (pi .* phi2 .* sqrt(phi2.^2 - phi1.^2));
			h_w = @(phi1,phi2)  min(h_Omega(phi1,phi2),h_Omega(phi1,phi1 + s_delta)).*(phi2<phi1+lambda_W/2);
			%%%
			v_clusterNum = 100:100:1000;

			ch_clustType = 'random';
			ch_estimationType = 'blind';
			s_kernelStd = 0.11;
			h_kernel = @(input1,input2) exp(-norm(input1-input2)./(2 * s_kernelStd^2));
			mu_w = 1e-5; % 1e-4 for l1, 1e-5 for total variation
			ini_F = rand(s_yAxiSize,s_xAxiSize);
						
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,~,~] = dataGenerator.realization();
			
			% B) estimation
			for s_itrIdx = 1:length(v_clusterNum)

			s_clusterNum = v_clusterNum(s_itrIdx);
			est = ChannelGainMapEstimator('mu_w',mu_w,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'ch_dataType',ch_dataType,'s_resolution',s_resolution,'s_clusterNum',s_clusterNum,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'lambda_W',lambda_W);
			[~,~,~,h_w_est] = est.optfunctionEstimation(t_sensorPos,t_sensorInd);

			% DISPLAY of weight functions
			s_lengthPhi1 = length(obj.v_Phi1real);
			s_lengthPhi2Axis = length(obj.v_Phi2Axisreal);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
			end
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1real,obj.v_rangPhi2real,obj.v_intervGridreal);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			fig_title=sprintf('Weight functions w/ N_c = %d',s_clusterNum);
			F(s_itrIdx) = F_figure('X',obj.v_Phi2Axisreal,'Y',m_evaluated_w,'colorp',12,'tit',fig_title,'leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			
			end
			
		end

		% This is a function to see changes in fitting on the normalized
		% ellipse weight function with the estimated weight function as the number
		% of cetroids increases.
		function F = compute_fig_6005(obj,niter)
						
			% Create data genator
			ch_dataType = 'real';
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.3937; % parameter to determine the threshold for nonzero weights
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			v_clusterNum = 500:100:1200;

			ch_clustType = 'random';
			ch_estimationType = 'blind';
			s_kernelStd = 0.11;
			h_kernel = @(input1,input2) exp(-norm(input1-input2).^2./(2 * s_kernelStd^2));
			mu_w = 1e-5; % 1e-4 for l1, 1e-5 for total variation
			ini_F = rand(s_yAxiSize,s_xAxiSize);
						
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,~,~] = dataGenerator.realization();
			
			% B) estimation
			for s_itrIdx = 1:length(v_clusterNum)

			s_clusterNum = v_clusterNum(s_itrIdx);
			est = ChannelGainMapEstimator('mu_w',mu_w,'h_w',h_w,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'ch_dataType',ch_dataType,'s_resolution',s_resolution,'s_clusterNum',s_clusterNum,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'lambda_W',lambda_W);
			[~,~,~,h_w_est] = est.optfunctionEstimation(t_sensorPos,t_sensorInd);

			% DISPLAY of weight functions
			s_lengthPhi1 = length(obj.v_Phi1real);
			s_lengthPhi2Axis = length(obj.v_Phi2Axisreal);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
			end
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1real,obj.v_rangPhi2real,obj.v_intervGridreal);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			fig_title=sprintf('Weight functions w/ N_c = %d',s_clusterNum);
			F(s_itrIdx) = F_figure('X',obj.v_Phi2Axisreal,'Y',m_evaluated_w,'colorp',12,'tit',fig_title,'leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			
			end
			
		end
		
		% This is a function to estimate optimal coefficient to
		% represent the weight function of the normalized ellipse model,
		% and interpolate the field by solving the original formulation in
		% the paper w.r.t the field, given the coefficients bar{alpha}.
		function F = compute_fig_6006(obj,niter)
						
			% Create data genator
			ch_dataType = 'real';
			dataGenerator = RealSensorMeasurementsGenerator();
			
			% Create estimator
			s_resolution = 1.5;
			s_widthY = 21;
			s_widthX = 21;
			s_xAxiSize = (s_widthX-1) * s_resolution + 1;
			s_yAxiSize = (s_widthY-1) * s_resolution + 1;
			lambda_W = 0.3937; % parameter to determine the threshold for nonzero weights
			h_w = @(phi1,phi2) (1/sqrt(phi1)).*(phi2<phi1+lambda_W/2); % normalized ellipse model
			s_clusterNum = 3000;
			ch_clustType = 'random';
			ch_estimationType = 'non-blind-non-para';
			ch_calibrationType = 'previous'; % 'none';'previous';'simultaneous'
			ch_reg_f_type = 'tikhonov'; %'l1_PCO'; %'totalvariation';
			s_kernelStd = 0.07;
			h_kernel = @(input1,input2) exp(-norm(input1-input2).^2./(2 * s_kernelStd^2));

			mu_f = 8 * 1e-4;
			mu_w = 1e-3; % 1e-4 for l1, 1e-5 for total variation
			rho =  1e-4; % for ISTA, rho is roughtly 2.5
			ini_F = rand(s_yAxiSize,s_xAxiSize);
			est = ChannelGainMapEstimator('mu_f',mu_f,'mu_w',mu_w,'h_w',h_w,'rho',rho,'ini_F',ini_F,'ch_estimationType',ch_estimationType,'ch_reg_f_type',ch_reg_f_type,'ch_dataType',ch_dataType,'ch_calibrationType',ch_calibrationType,'s_resolution',s_resolution,'s_clusterNum',s_clusterNum,'ch_clustType',ch_clustType,'h_kernel',h_kernel,'lambda_W',lambda_W);
						
			% SIMULATION
			% A) data generation
			[t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = dataGenerator.realization();
			
			% B) estimation
			[m_F_est,h_w_est] = est.estimate(t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing);

			% DISPLAY of weight functions
			s_lengthPhi1 = length(obj.v_Phi1real);
			s_lengthPhi2Axis = length(obj.v_Phi2Axisreal);
			for s_phi1Ind = 1 : s_lengthPhi1
				leg{s_phi1Ind} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
				leg{s_phi1Ind + s_lengthPhi1} = sprintf('\\phi_1 = %g', obj.v_Phi1real(s_phi1Ind));
			end
			[m_w_o,m_w_hat] = ChannelGainMapEstimator.evaluate_w(h_w,h_w_est,obj.v_rangPhi1real,obj.v_rangPhi2real,obj.v_intervGridreal);
			m_evaluated_w = [m_w_o(1:12,1:s_lengthPhi2Axis);m_w_hat(1:12,1:s_lengthPhi2Axis)];
			F(1) = F_figure('X',obj.v_Phi2Axisreal,'Y',m_evaluated_w,'colorp',12,'tit','Weight functions','leg',leg,'xlab','\phi_2','ylab','FUNCTION VALUE');
			F(2) = F_figure('Z',ChannelGainMapEstimator.postprocessReal(m_F_est),'tit','Estimated','pos',[100.0000  402.0000  [350.3077  350.3077]]);

		end
		
	end
	
	
end
