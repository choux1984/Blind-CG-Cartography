classdef ChannelGainMapEstimator < Parameter

	properties(Constant)
		eps_error = 1e-10;
		max_itr = 200;
	end
		
	properties % Required by superclass Parameter
		c_parsToPrint    = {};
		c_stringToPrint  = {};
		c_patternToPrint = {};
	end
	
	properties
		ini_F        % Nx x Ny matrix with the initialization for the spatial loss field
		             % Mandatory in all cases since this is what decides
		             % the size of the estimate of F
		ch_reg_f_type% regularization type for the SLF:
		             %     'tikhonov'= Tikhonov;
		             %     'l1_ISTA' = l_1 solved with ISTA; 'l1_PCO' = l_1 solved with pathwise coordinate opt.;
		             %     'totalvariation' = total variation
		mu_f         % Regularization parameter for spatial loss field
				
		ch_estimationType
		             % 'non-blind'
					 % 'blind'
					 % 'warm'
					 % 'blind-only-ellipse'
		
        % calibration
		ch_calibrationType = 'none';
		             %  'none'  : values of pathloss and gains taken from
		             %            the properties below
					 %  'simultaneous' 
		v_gains;     % (n_sensor_types x 1) vector, where the i-th entry 
		             % contains the gain of the i-th type of sensor. 
		s_pathLossExponent;
		
		m_Omega = []; % (n_measurements x n_sensor_types) matrix with a 1
		             % in the (n,m) position if a sensor of type m is the tx
		             % or rx of measurement n
					 %
					 % by default it is assumed that all sensors are of the
					 % same type, so m_Omega = ones(n_sensor_types,1).
		
		
		%ch_gainType = 'different';
		             % parameter used when ch_calibrationType =
		             % 'simultaneous'
		             %  'different' : a gain value is estimated per sensor
					 %  'same'      : a gain value is estimated for all
					 %                sensors
					 % 'custom' : one gain value is estimated per column of 
					 %            matrix m_Omega.
					 
	    			 
		% estimation
		rho          % step size for (ISTA/ADMM) algorithm
		s_resolution = 1;% scaling factor to decide a resolution of an output image, e.g s_resolution = 1 (default)
		                % resolution of the output changes only for real dataset.    
		
		 
		% non-blind estimation
		h_w = [];    %  
		
		% blind estimation
		mu_w         % Regularization parameter for weight function		
		h_kernel;    % Kernel function to interpolate a weight function. e.g. Gaussian kernel
		s_clusterNum % Number of clusters; Set s_clusterNum = [] for no clustering
		ch_clustType % Types of feature clusterting, and corresponding cluster centroids. 
		             % 'random' : random sampling 
					 % 'kmeans' : k-means
		lambda_W     % h_w is estimated with domain equal to an ellipse whose semi-minor axis is lambda_W
		ch_coefficient % name of a file containing kernel coefficients and centroids used to obtain it
		s_SemiAxisLength4Sample % Parameter to define the support for feature-sampling to find centroids.
	end
		
	methods % estimator
		
		function obj = ChannelGainMapEstimator(varargin)
			obj@Parameter(varargin{:});
		end
				
		function [m_F_est,h_w_est] = estimate(obj,m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing)	
			%
			% OUTPUT:
			%   m_est_F         N_x-by-N_y matrix with the estimate of F
			%   h_w_est    - non-blind estimation
			%                 h_w_est = []
			%              - blind estimation
			%                 h_w_est is the estimate of h_w_est
			%
			%  m_centroids     2-by-(n_measurements*Ng): collection of all feature vectors
			
			
			% dependent variables
			[N_x, N_y] = size(obj.ini_F);
						
			% Estimation
			switch obj.ch_estimationType
				case 'non-blind'					
					m_F_est  = obj.nonblindEstimation(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);										
					h_w_est = [];
				case 'non-blind-non-para'
					[m_F_est,h_w_est]  = obj.nonblindEstimationNonParametric(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
				case 'blind'
					[m_F_est,h_w_est] = obj.blindEstimation(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);			
				case 'warm'
					[m_F_est,h_w_est] = obj.blindEstimationWarmStart(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
				case 'blind-only-ellipse'
					[m_F_est,h_w_est] = obj.blindEstimationOnlyEllipse(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);			
				otherwise
					error('unrecognized option');
			end
				 
		end
		
		function m_F_est  = nonblindEstimation(obj,t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing)
			% Estimate a spatial loss field in non-blind fashion with a
			% chosen regularizer (tikhonov, l-1, and total variation).
			
			% With slight abuse of notation, the argument m_sensorPos
			% should be a tensor for "real" dataset (1st slab for the
			% non-structure test and 2nd slab for the stuctured test since
			% locations of sensors in both tests are not coincident due to
			% the mislaignment of sensor locations). 
								         
			% check
			if isempty(obj.ini_F)&&(strcmp(obj.ch_reg_f_type,'tikhonov')==0)
				assert(~isempty(obj.ini_F));				
			end
			
			if strcmp(obj.ch_dataType,'real')
				assert(~strcmp(obj.ch_calibrationType,'none'),'Selected option is not valid.');
			end
            
            % 1. Compute a weight matrix for each pair of Tx / Rx from m_txPos / m_rxPos
            [N_x, N_y] = size(obj.ini_F);
			s_gridNum = N_x * N_y;
            s_measurementNum = length(v_measurements); 
			x_axis = repmat((0:N_y-1)./obj.s_resolution, N_x, 1);   % x_axis of a grid
			y_axis = repmat((N_x-1:-1:0)'./obj.s_resolution, 1, N_y); % y_axis of a grid
			
			switch obj.ch_dataType
				case 'synthetic'
					m_sensorPos=t_sensorPos;
					m_sensorInd=t_sensorInd;
				case 'real'
					m_sensorPos=t_sensorPos(:,:,2);
					m_sensorInd=t_sensorInd(:,:,2);
					s_measurementNum = s_measurementNum + 20; % compensate 20 missing measurements
			end
			
            m_vecWCollection = zeros(s_gridNum,s_measurementNum);
			for s_measurementInd = 1 : s_measurementNum
				v_txPos = m_sensorPos(:,m_sensorInd(1,s_measurementInd));
				v_rxPos = m_sensorPos(:,m_sensorInd(2,s_measurementInd));
				
				s_phi1 = norm(v_txPos-v_rxPos);
				m_phi2 = sqrt( (x_axis-v_txPos(1)).^2 + (y_axis-v_txPos(2)).^2 ) + sqrt( (x_axis-v_rxPos(1)).^2 + (y_axis-v_rxPos(2)).^2 );
				
% 				m_W = obj.h_w(s_phi1,m_phi2);
				m_W = zeros(size(m_phi2,1),size(m_phi2,2));
				for s_ygridInd = 1 : size(m_phi2,1)
					for s_xgridInd = 1 : size(m_phi2,2)
						if m_phi2(s_ygridInd,s_xgridInd) < s_phi1 + 0.5 * obj.lambda_W
							m_W(s_ygridInd,s_xgridInd) = obj.h_w(s_phi1,m_phi2(s_ygridInd,s_xgridInd));
						end
					end
				end
				
				m_vecWCollection(:,s_measurementInd) = m_W(:);
			end
			
            % 2. Estimate m_F_est according to obj.ch_reg_f_type (regularizer type) and ch_calibrationType.
			switch obj.ch_calibrationType
				case 'none'
					[v_sumOfGains,v_pathLoss] = obj.estimateSensorGainAndPathLoss(m_sensorPos,m_sensorInd,[]);
					s_check = -1 * (v_measurements - (v_sumOfGains - v_pathLoss));
					v_f_est = obj.chooseSolver(s_check,m_vecWCollection);
				case 'previous'
					% for real data simulation, m_sensorPos' for structured
					% and free spaces are required
					switch obj.ch_dataType
						case 'synthetic'
							[v_sumOfGains,v_pathLoss] = obj.estimateSensorGainAndPathLoss(m_sensorPos,m_sensorInd,v_measurementsNoShadowing);
						case 'real'
							[v_sumOfGains,v_pathLoss] = obj.parameterEstimatorReal(t_sensorPos,t_sensorInd,v_measurementsNoShadowing);
							m_vecWCollection = m_vecWCollection(:,[1:580,601:2400]); % compensate 20 missing measurements
					end
					s_check =  -1 * (v_measurements - (v_sumOfGains - v_pathLoss));
					v_f_est = obj.chooseSolver(s_check,m_vecWCollection);
				case 'simultaneous'
					% we need an estimator jointly estimating v_f_est,
					% pathloss exponent, and tx/rx gains only with shadow-faded channel
					% gain measurements

					v_sensorDistances = zeros(s_measurementNum,1);
					for s_measurementInd = 1: s_measurementNum
						v_txPos = m_sensorPos(:,m_sensorInd(1,s_measurementInd));
						v_rxPos = m_sensorPos(:,m_sensorInd(2,s_measurementInd));
						v_sensorDistances(s_measurementInd) = 10 * log10(norm(v_txPos-v_rxPos));
					end
					
					switch obj.ch_dataType
						case 'synthetic'
							s_sensorNum = size(m_sensorPos,2);
							m_Omega = obj.sensorMapOp(m_sensorInd,s_sensorNum);
						case 'real'
							s_sensorNum = 120;
							s_measurementNum = 2380;
							m_Omega = csvread('m_Omega.csv'); % due to measurements from sensors not conforming to the sensor selction rule
							m_Omega = m_Omega([1:580,601:2400],:); % remove rows corresponding to messing measurements
							m_vecWCollection = m_vecWCollection(:,[1:580,601:2400]); % remove cols corresponding to messing measurements
							v_sensorDistances = v_sensorDistances([1:580,601:2400],:);
					end
									
					m_barI = eye(s_measurementNum) - (1/(norm(v_sensorDistances,2)^2))*(v_sensorDistances)*(v_sensorDistances');
					m_tempOmega = m_barI * m_Omega;
					% add "1e-4 * eye(s_sensorNum)" to "m_tempOmega' * m_tempOmega" for matrix inversion
					m_barOmega = m_Omega * ((m_tempOmega' * m_tempOmega + 1e-4 * eye(s_sensorNum)) \ (m_tempOmega' * m_barI));
					m_tildeI = m_barI * (eye(s_measurementNum) - m_barOmega);
					v_barmeasurements = m_tildeI * v_measurements;
					m_barvecWCollection = -1 * m_vecWCollection * m_tildeI';
					v_f_est = obj.chooseSolver(v_barmeasurements,m_barvecWCollection);
			end
          
			m_F_est = reshape(v_f_est,N_x,N_y);			
		end
		
		function [m_F_est,h_w_est]  = nonblindEstimationNonParametric(obj,t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing)
			% Estimate a spatial loss field in non-blind fashion with a
			% chosen regularizer (tikhonov, l-1, and total variation).
			% Different to the nonblind estimation, m_F_est is estimated by
			% using vectorized weight matrices represented by the kernel and
			% mapping matrices with v_alpha optimally obatined from solving a kernel
			% regression problem (see m_A in the paper).
			
			% With slight abuse of notation, the argument m_sensorPos
			% should be a tensor for "real" dataset (1st slab for the
			% non-structure test and 2nd slab for the stuctured test since
			% locations of sensors in both tests are not coincident due to
			% the mislaignment of sensor locations). 
								         
			% check
			if isempty(obj.ini_F)&&(strcmp(obj.ch_reg_f_type,'tikhonov')==0)
				assert(~isempty(obj.ini_F));				
			end
			
			if strcmp(obj.ch_dataType,'real')
				assert(~strcmp(obj.ch_calibrationType,'none'),'Selected option is not valid.');
			end
            
            % 1. Compute a weight matrix for each pair of Tx / Rx from m_txPos / m_rxPos
            [N_x, N_y] = size(obj.ini_F);
			s_gridNum = N_x * N_y;
            s_measurementNum = length(v_measurements); 
			x_axis = repmat((0:N_y-1)./obj.s_resolution, N_x, 1);   % x_axis of a grid
			y_axis = repmat((N_x-1:-1:0)'./obj.s_resolution, 1, N_y); % y_axis of a grid
			
			switch obj.ch_dataType
				case 'synthetic'
					m_sensorPos=t_sensorPos;
					m_sensorInd=t_sensorInd;
				case 'real'
					m_sensorPos=t_sensorPos(:,:,2);
					m_sensorInd=t_sensorInd(:,:,2);
					s_measurementNum = s_measurementNum + 20; % compensate 20 missing measurements
			end
			
			[m_centroids, v_centroidsInd,v_alpha,h_w_est] = obj.optfunctionEstimation(t_sensorPos,t_sensorInd);
			m_K = ChannelGainMapEstimator.kernelMatrix(m_centroids,obj.h_kernel);
			m_RK = m_K(v_centroidsInd,:); % precomputation of "R*m_K"
			v_RKa = m_RK * v_alpha;
			m_vecWCollection = reshape(v_RKa,s_gridNum,s_measurementNum); % see the definition of "A" in the paper

			
            % 2. Estimate m_F_est according to obj.ch_reg_f_type (regularizer type) and ch_calibrationType.
			switch obj.ch_calibrationType
				case 'none'
					[v_sumOfGains,v_pathLoss] = obj.estimateSensorGainAndPathLoss(m_sensorPos,m_sensorInd,[]);
					s_check = -1 * (v_measurements - (v_sumOfGains - v_pathLoss));
					v_f_est = obj.chooseSolver(s_check,m_vecWCollection);
				case 'previous'
					% for real data simulation, m_sensorPos' for structured
					% and free spaces are required
					switch obj.ch_dataType
						case 'synthetic'
							[v_sumOfGains,v_pathLoss] = obj.estimateSensorGainAndPathLoss(m_sensorPos,m_sensorInd,v_measurementsNoShadowing);
						case 'real'
							[v_sumOfGains,v_pathLoss] = obj.parameterEstimatorReal(t_sensorPos,t_sensorInd,v_measurementsNoShadowing);
							m_vecWCollection = m_vecWCollection(:,[1:580,601:2400]); % compensate 20 missing measurements
					end
					s_check =  -1 * (v_measurements - (v_sumOfGains - v_pathLoss));
					v_f_est = obj.chooseSolver(s_check,m_vecWCollection);
				case 'simultaneous'
					% we need an estimator jointly estimating v_f_est,
					% pathloss exponent, and tx/rx gains only with shadow-faded channel
					% gain measurements

					v_sensorDistances = zeros(s_measurementNum,1);
					for s_measurementInd = 1: s_measurementNum
						v_txPos = m_sensorPos(:,m_sensorInd(1,s_measurementInd));
						v_rxPos = m_sensorPos(:,m_sensorInd(2,s_measurementInd));
						v_sensorDistances(s_measurementInd) = 10 * log10(norm(v_txPos-v_rxPos));
					end
					
					switch obj.ch_dataType
						case 'synthetic'
							s_sensorNum = size(m_sensorPos,2);
							m_Omega = obj.sensorMapOp(m_sensorInd,s_sensorNum);
						case 'real'
							s_sensorNum = 120;
							s_measurementNum = 2380;
							m_Omega = csvread('m_Omega.csv'); % due to measurements from sensors not conforming to the sensor selction rule
							m_Omega = m_Omega([1:580,601:2400],:); % remove rows corresponding to messing measurements
							m_vecWCollection = m_vecWCollection(:,[1:580,601:2400]); % remove cols corresponding to messing measurements
							v_sensorDistances = v_sensorDistances([1:580,601:2400],:);
					end
									
					m_barI = eye(s_measurementNum) - (1/(norm(v_sensorDistances,2)^2))*(v_sensorDistances)*(v_sensorDistances');
					m_tempOmega = m_barI * m_Omega;
					% add "1e-4 * eye(s_sensorNum)" to "m_tempOmega' * m_tempOmega" for matrix inversion
					m_barOmega = m_Omega * ((m_tempOmega' * m_tempOmega + 1e-4 * eye(s_sensorNum)) \ (m_tempOmega' * m_barI));
					m_tildeI = m_barI * (eye(s_measurementNum) - m_barOmega);
					v_barmeasurements = m_tildeI * v_measurements;
					m_barvecWCollection = -1 * m_vecWCollection * m_tildeI';
					v_f_est = obj.chooseSolver(v_barmeasurements,m_barvecWCollection);
			end
          
			m_F_est = reshape(v_f_est,N_x,N_y);			
		end
		
		function [m_F_est,h_w_est] = blindEstimation(obj,t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing)
			% Estimate a spatial loss field in blind fashion with a
			% chosen regularizer for the SLF(tikhonov, l-1, and total variation).
			% In addition, a weight function is estimated with a kernel
			% method. In particular, the weight function lies in
			% reproducing kernel hilbert space. Complexity of the function
			% is controllable by a hilbert norm.
            
			if strcmp(obj.ch_dataType,'real')
				assert(~strcmp(obj.ch_calibrationType,'none'),'Selected option is not valid.');
			end
						
            % 1. Compute features and correponding kernel matrix. For a
            % practical purpose, features are clustered and represented by
            % much less number of cluster centroids. Those
            % centroids are found by either 'random' sampling, or
            % well-known 'k-means'.
			
            [N_x, N_y] = size(obj.ini_F);
            s_gridNum = N_x * N_y;
            s_measurementNum = length(v_measurements);

			switch obj.ch_dataType
				case 'synthetic'
					m_sensorPos=t_sensorPos;
					m_sensorInd=t_sensorInd;
				case 'real'
					m_sensorPos=t_sensorPos(:,:,2);
					m_sensorInd=t_sensorInd(:,:,2);
					s_measurementNum = s_measurementNum + 20; % compensate 20 missing measurements
			end
			
			[m_centroids, v_centroidsInd]= obj.findCentroids(m_sensorPos,m_sensorInd);
			m_K = ChannelGainMapEstimator.kernelMatrix(m_centroids,obj.h_kernel);
			
            % 2. Estimate m_F_est according to obj.ch_reg_f_type (regularizer type) and ch_calibrationType
			%    and v_alpha.
				
			% alternating minimization setup
            s_stoppingCriterion = 1e-5;
			s_iterationMax = 1e2;
            estimateDifference = inf;
			m_RK = m_K(v_centroidsInd,:); % precomputation of "R*m_K"
			prev_v_f = obj.ini_F(:);

            switch obj.ch_calibrationType
				case 'none'
					[v_sumOfGains,v_pathLoss] = obj.estimateSensorGainAndPathLoss(m_sensorPos,m_sensorInd,[]);
					s_check = -1 * (v_measurements - (v_sumOfGains - v_pathLoss));
					
					s_iterationNum = 1;
					while (estimateDifference > s_stoppingCriterion) && (s_iterationNum < s_iterationMax)
						
						% [S1] estimate v_alpha
						m_IfR = obj.computeIfR(prev_v_f,v_centroidsInd);
						m_IfRK = m_IfR * m_K;										
						v_alpha = ChannelGainMapEstimator.ridgeRegression(m_IfRK',s_check,obj.mu_w,m_K);
						
						% [S2] estimate v_f_est
						v_RKa = m_RK * v_alpha;
						m_A = reshape(v_RKa,s_gridNum,s_measurementNum)'; % see the definition of "A" in the paper
						v_f_est = obj.chooseSolver(s_check,m_A');
	
						estimateDifference = norm(prev_v_f - v_f_est,2);
						prev_v_f = v_f_est;
						s_iterationNum = s_iterationNum + 1;
					end
								
				case 'previous'
					switch obj.ch_dataType
						case 'synthetic'
							[v_sumOfGains,v_pathLoss] = obj.estimateSensorGainAndPathLoss(m_sensorPos,m_sensorInd,v_measurementsNoShadowing);
							s_check = -1 * (v_measurements - (v_sumOfGains - v_pathLoss));
						case 'real'
							[v_sumOfGains,v_pathLoss] = obj.parameterEstimatorReal(t_sensorPos,t_sensorInd,v_measurementsNoShadowing);
							s_check = -1 * (v_measurements - (v_sumOfGains - v_pathLoss));
% 							s_check = [s_check(1:580,1);zeros(20,1);s_check(581:end,1)];
					end
					s_iterationNum = 1;

					while (estimateDifference > s_stoppingCriterion) && (s_iterationNum < s_iterationMax)			
						% [S1] estimate v_alpha
						m_IfR = obj.computeIfR(prev_v_f,v_centroidsInd);
						m_IfRK = m_IfR * m_K;		
						if strcmp(obj.ch_dataType,'real')
							m_IfRK = m_IfRK([1:580,601:2400],:);
						end
						
						v_alpha = ChannelGainMapEstimator.ridgeRegression(m_IfRK',s_check,obj.mu_w,m_K);
						
						% [S2] estimate v_f_est
						v_RKa = m_RK * v_alpha;
						m_A = reshape(v_RKa,s_gridNum,s_measurementNum)'; % see the definition of "A" in the paper						
						if strcmp(obj.ch_dataType,'real')
							m_A = m_A([1:580,601:2400],:);
						end
						
						v_f_est = obj.chooseSolver(s_check,m_A');
	
						estimateDifference = norm(prev_v_f - v_f_est,2);
						prev_v_f = v_f_est;
						s_iterationNum = s_iterationNum + 1;
					end
				case 'simultaneous'
					% we need an estimator jointly estimating v_f_est,
					% pathloss exponent, and tx/rx gains only with shadow-faded channel
					% gain measurements
					s_sensorNum = size(m_sensorPos,2);
					v_sensorDistancesdB = zeros(s_measurementNum,1);
					
					if strcmp(obj.ch_gainType,'different') == 1
						m_Omega = obj.sensorMapOp(m_sensorInd,s_sensorNum);
					else
						m_Omega = ones(s_measurementNum,1);
					end

					for s_measurementInd = 1: s_measurementNum
						v_txPos = m_sensorPos(:,m_sensorInd(1,s_measurementInd));
						v_rxPos = m_sensorPos(:,m_sensorInd(2,s_measurementInd));
						v_sensorDistancesdB(s_measurementInd) = 10 * log10(norm(v_txPos-v_rxPos));
					end
					
					m_barI = eye(s_measurementNum) - (1/(norm(v_sensorDistancesdB,2)^2))*(v_sensorDistancesdB)*(v_sensorDistancesdB');
					m_tempOmega = m_barI * m_Omega;
					% add "1e-4 * eye(s_sensorNum)" to "m_tempOmega' * m_tempOmega" for matrix inversion
					m_barOmega = m_Omega * ((m_tempOmega' * m_tempOmega + 1e-4 * eye(size(m_Omega,2))) \ (m_tempOmega' * m_barI));
					m_tildeI = m_barI * (eye(s_measurementNum) - m_barOmega );
					v_barmeasurements = m_tildeI * v_measurements;
										
					s_iterationNum = 1;
					while (estimateDifference > s_stoppingCriterion) && (s_iterationNum < s_iterationMax)
						% [S1] estimate v_alpha
						m_IfR = obj.computeIfR(prev_v_f,v_centroidsInd);
						m_IfRK = m_IfR * m_K;
						m_barIfKR = -1 * m_IfRK' * m_tildeI';
						v_alpha = ChannelGainMapEstimator.ridgeRegression(m_barIfKR,v_barmeasurements,obj.mu_w,m_K);
						
						% [S2] estimate v_f_est
						v_RKa = m_RK * v_alpha;
						m_A = reshape(v_RKa,s_gridNum,s_measurementNum)'; % see the definition of "A" in the paper
						m_barA = -1 * m_A' * m_tildeI';
						v_f_est = obj.chooseSolver(v_barmeasurements,m_barA);
						
						estimateDifference = norm(prev_v_f - v_f_est,2);
						prev_v_f = v_f_est;
						s_iterationNum = s_iterationNum + 1;
					end
					
			end
          
			% blind estimator outputs
			m_F_est = reshape(v_f_est,N_x,N_y);						
			h_w_est = @(v_input1,v_input2) obj.represeThm(v_alpha,m_centroids,[v_input1;v_input2]);
			
		end
		
		function [m_F_est,h_w_est] = blindEstimationWarmStart(obj,t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing)
			% Estimate a spatial loss field in blind fashion with a
			% chosen regularizer for the SLF(tikhonov, l-1, and total variation).
			% In addition, a weight function is estimated with a kernel
			% method. In particular, the weight function lies in
			% reproducing kernel hilbert space. Complexity of the function
			% is controllable by a hilbert norm.
            
			if strcmp(obj.ch_dataType,'real')
				assert(~strcmp(obj.ch_calibrationType,'none'),'Selected option is not valid.');
			end
						
            % 1. Compute features and correponding kernel matrix. For a
            % practical purpose, features are clustered and represented by
            % much less number of cluster centroids. Those
            % centroids are found by either 'random' sampling, or
            % well-known 'k-means'.
			
            [N_x, N_y] = size(obj.ini_F);
            s_gridNum = N_x * N_y;
            s_measurementNum = length(v_measurements);

			switch obj.ch_dataType
				case 'synthetic'
					m_sensorPos=t_sensorPos;
					m_sensorInd=t_sensorInd;
				case 'real'
					m_sensorPos=t_sensorPos(:,:,2);
					m_sensorInd=t_sensorInd(:,:,2);
					s_measurementNum = s_measurementNum + 20; % compensate 20 missing measurements
			end
			
			load(obj.ch_coefficient);
			m_K = ChannelGainMapEstimator.kernelMatrix(m_centroids,obj.h_kernel);
			
            % 2. Estimate m_F_est according to obj.ch_reg_f_type (regularizer type) and ch_calibrationType
			%    and v_alpha.
				
			% alternating minimization setup
            s_stoppingCriterion = 1e-5;
			s_iterationMax = 2;
            estimateDifference = inf;
			m_RK = m_K(v_centroidsInd,:); % precomputation of "R*m_K"
			
			switch obj.ch_calibrationType
				case 'none'
					[v_sumOfGains,v_pathLoss] = obj.estimateSensorGainAndPathLoss(m_sensorPos,m_sensorInd,[]);
					s_check = -1 * (v_measurements - (v_sumOfGains - v_pathLoss));
					
					s_iterationNum = 1;
					while (estimateDifference > s_stoppingCriterion) && (s_iterationNum < s_iterationMax)
						
						% [S1] estimate v_alpha
						m_IfR = obj.computeIfR(prev_v_f,v_centroidsInd);
						m_IfRK = m_IfR * m_K;
						v_alpha = ChannelGainMapEstimator.ridgeRegression(m_IfRK',s_check,obj.mu_w,m_K);
						
						% [S2] estimate v_f_est
						v_RKa = m_RK * v_alpha;
						m_A = reshape(v_RKa,s_gridNum,s_measurementNum)'; % see the definition of "A" in the paper
						v_f_est = obj.chooseSolver(s_check,m_A');
						
						estimateDifference = norm(prev_v_f - v_f_est,2);
						prev_v_f = v_f_est;
						s_iterationNum = s_iterationNum + 1;
					end
					
				case 'previous'
					switch obj.ch_dataType
						case 'synthetic'
							[v_sumOfGains,v_pathLoss] = obj.estimateSensorGainAndPathLoss(m_sensorPos,m_sensorInd,v_measurementsNoShadowing);
							s_check = -1 * (v_measurements - (v_sumOfGains - v_pathLoss));
						case 'real'
							[v_sumOfGains,v_pathLoss] = obj.parameterEstimatorReal(t_sensorPos,t_sensorInd,v_measurementsNoShadowing);
							s_check = -1 * (v_measurements - (v_sumOfGains - v_pathLoss));
					end
					s_iterationNum = 1;
					
					% Test
					v_RKa = m_RK * v_alpha;
					m_A = reshape(v_RKa,s_gridNum,s_measurementNum)'; % see the definition of "A" in the paper
					if strcmp(obj.ch_dataType,'real')
						m_A = m_A([1:580,601:2400],:);
					end
					
					prev_v_f = obj.chooseSolver(s_check,m_A');
					
					
					while (estimateDifference > s_stoppingCriterion) && (s_iterationNum < s_iterationMax)
						% [S1] estimate v_alpha
						m_IfR = obj.computeIfR(prev_v_f,v_centroidsInd);
						m_IfRK = m_IfR * m_K;
						if strcmp(obj.ch_dataType,'real')
							m_IfRK = m_IfRK([1:580,601:2400],:);
						end
						
						v_alpha = ChannelGainMapEstimator.ridgeRegression(m_IfRK',s_check,obj.mu_w,m_K);
						
						% [S2] estimate v_f_est
						v_RKa = m_RK * v_alpha;
						m_A = reshape(v_RKa,s_gridNum,s_measurementNum)'; % see the definition of "A" in the paper
						if strcmp(obj.ch_dataType,'real')
							m_A = m_A([1:580,601:2400],:);
						end
						
						v_f_est = obj.chooseSolver(s_check,m_A');
						
						estimateDifference = norm(prev_v_f - v_f_est,2);
						prev_v_f = v_f_est;
						s_iterationNum = s_iterationNum + 1;
					end
				case 'simultaneous'
					% we need an estimator jointly estimating v_f_est,
					% pathloss exponent, and tx/rx gains only with shadow-faded channel
					% gain measurements
					s_sensorNum = size(m_sensorPos,2);
					v_sensorDistancesdB = zeros(s_measurementNum,1);
					
					if strcmp(obj.ch_gainType,'different') == 1
						m_Omega = obj.sensorMapOp(m_sensorInd,s_sensorNum);
					else
						m_Omega = ones(s_measurementNum,1);
					end
					
					for s_measurementInd = 1: s_measurementNum
						v_txPos = m_sensorPos(:,m_sensorInd(1,s_measurementInd));
						v_rxPos = m_sensorPos(:,m_sensorInd(2,s_measurementInd));
						v_sensorDistancesdB(s_measurementInd) = 10 * log10(norm(v_txPos-v_rxPos));
					end
					
					m_barI = eye(s_measurementNum) - (1/(norm(v_sensorDistancesdB,2)^2))*(v_sensorDistancesdB)*(v_sensorDistancesdB');
					m_tempOmega = m_barI * m_Omega;
					% add "1e-4 * eye(s_sensorNum)" to "m_tempOmega' * m_tempOmega" for matrix inversion
					m_barOmega = m_Omega * ((m_tempOmega' * m_tempOmega + 1e-4 * eye(size(m_Omega,2))) \ (m_tempOmega' * m_barI));
					m_tildeI = m_barI * (eye(s_measurementNum) - m_barOmega );
					v_barmeasurements = m_tildeI * v_measurements;
					
					s_iterationNum = 1;
					while (estimateDifference > s_stoppingCriterion) && (s_iterationNum < s_iterationMax)
						% [S1] estimate v_alpha
						m_IfR = obj.computeIfR(prev_v_f,v_centroidsInd);
						m_IfRK = m_IfR * m_K;
						m_barIfKR = -1 * m_IfRK' * m_tildeI';
						v_alpha = ChannelGainMapEstimator.ridgeRegression(m_barIfKR,v_barmeasurements,obj.mu_w,m_K);
						
						% [S2] estimate v_f_est
						v_RKa = m_RK * v_alpha;
						m_A = reshape(v_RKa,s_gridNum,s_measurementNum)'; % see the definition of "A" in the paper
						m_barA = -1 * m_A' * m_tildeI';
						v_f_est = obj.chooseSolver(v_barmeasurements,m_barA);
						
						estimateDifference = norm(prev_v_f - v_f_est,2);
						prev_v_f = v_f_est;
						s_iterationNum = s_iterationNum + 1;
					end
					
			end
			
			% blind estimator outputs
			m_F_est = reshape(v_f_est,N_x,N_y);
			h_w_est = @(v_input1,v_input2) obj.represeThm(v_alpha,m_centroids,[v_input1;v_input2]);
			
		end

		function [m_F_est,h_w_est] = blindEstimationOnlyEllipse(obj,t_sensorPos,t_sensorInd,v_measurements)
			% Estimate a spatial loss field in blind fashion with a
			% chosen regularizer for the SLF(tikhonov, l-1, and total variation).
			% In addition, a weight function is estimated with a kernel
			% method. In particular, the weight function lies in
			% reproducing kernel hilbert space. Complexity of the function
			% is controllable by a hilbert norm. 
			
			% Furthermore, in this simulation, we only consider spatial
			% loss field entries only within an ellipse (accordingly other
			% related parameters, feature vectors, etc.)
			
% 			if strcmp(obj.ch_dataType,'real')
% 				assert(~strcmp(obj.ch_calibrationType,'none'),'Selected option is not valid.');
% 			end
			
			% 1. Compute features and correponding kernel matrix. For a
			% practical purpose, features are clustered and represented by
			% much fewer cluster centroids. Those
			% centroids are found by either 'random' sampling, or
			% well-known 'k-means'.
			
			[N_x, N_y] = size(obj.ini_F);
			s_gridNum = N_x * N_y;
			s_measurementNum = length(v_measurements);
			
% 			switch obj.ch_dataType
% 				case 'synthetic'
% 					m_sensorPos=t_sensorPos;
% 					m_sensorInd=t_sensorInd;
% 				case 'real'
% 					m_sensorPos=t_sensorPos(:,:,2);
% 					m_sensorInd=t_sensorInd(:,:,2);
% 					s_measurementNum = s_measurementNum + 20; % compensate 20 missing measurements
% 			end
			
			% It is required to find vectors cotaining indices of nonzero w (for example, 1 for the 
			% first element of a vectorized W matrix, and L for the last element), a size
			% of nonzero weights for each pair of sensors.
			[m_centroids, v_centroidsInd, v_nzWIdx,v_nzWbinSize]= obj.findCentroidsEllip(m_sensorPos,m_sensorInd);
			v_nz_PsiRKaIdx = ChannelGainMapEstimator.nzRelativeIdxFinder(v_nzWIdx,v_nzWbinSize,s_gridNum);
			m_K = ChannelGainMapEstimator.kernelMatrix(m_centroids,obj.h_kernel);
			
			% 2. Estimate m_F_est according to obj.ch_reg_f_type (regularizer type) and ch_calibrationType
			%    and v_alpha. As we only consider the grids corresponding
			%    to nonzero weights, m_IfRK and m_A should be changed
			%    by using v_nzWIdx and v_nzWbinSize.
			
	v_sensorDistances = [];
	        if isempty(obj.m_Omega)
				obj.m_Omega = ones(s_measurementNum,1);
				assert(length(obj.v_gains)==1);
			end
			
			% alternating minimization setup
			s_stoppingCriterion = 1e-5;
			s_iterationMax = 2 * 1e2;
			estimateDifference = inf;
			m_RK = m_K(v_centroidsInd,:); % precomputation of "R*m_K"
			prev_v_f = obj.ini_F(:);
			
			switch obj.ch_calibrationType
				case 'none'
					%[v_sumOfGains,v_pathLoss] = obj.estimateSensorGainAndPathLoss(m_sensorPos,m_sensorInd,[]);
					s_check = -1 * (v_measurements - (obj.m_Omega*obj.v_gains - obj.s_pathLossExponent*10*log10(v_sensorDistances)));
					
					s_iterationNum = 1;
					while (estimateDifference > s_stoppingCriterion) && (s_iterationNum < s_iterationMax)
						
						% [S1] estimate v_alpha
						m_IfPsiRK = obj.computeIfPsiRK(prev_v_f,v_nzWIdx,v_nzWbinSize,m_RK);			
						v_alpha = ChannelGainMapEstimator.ridgeRegression(m_IfPsiRK',s_check,obj.mu_w,m_K);
						
						% [S2] estimate v_f_est
						v_RKa = m_RK * v_alpha;
						v_tempA = zeros(s_gridNum * s_measurementNum,1);
						v_tempA(v_nz_PsiRKaIdx) = v_RKa;
						m_A = reshape(v_tempA,s_gridNum,s_measurementNum)'; % see the definition of "A" in the paper
						v_f_est = obj.chooseSolver(s_check,m_A');
						
						estimateDifference = norm(prev_v_f - v_f_est,2);
						prev_v_f = v_f_est;
						s_iterationNum = s_iterationNum + 1;
					end
					
% 				case 'previous'
% 					switch obj.ch_dataType
% 						case 'synthetic'
% 							[v_sumOfGains,v_pathLoss] = obj.estimateSensorGainAndPathLoss(m_sensorPos,m_sensorInd,v_measurementsNoShadowing);
% 							s_check = -1 * (v_measurements - (v_sumOfGains - v_pathLoss));
% 						case 'real'
% 							[v_sumOfGains,v_pathLoss] = obj.estimateSensorGainAndPathLossReal(t_sensorPos,t_sensorInd,v_measurementsNoShadowing);
% 							s_check = -1 * (v_measurements - (v_sumOfGains - v_pathLoss));
% 							% 							s_check = [s_check(1:580,1);zeros(20,1);s_check(581:end,1)];
% 					end
% 					s_iterationNum = 1;
% 					
% 					while (estimateDifference > s_stoppingCriterion) && (s_iterationNum < s_iterationMax)
% 						% [S1] estimate v_alpha
% 						m_IfPsiRK = obj.computeIfPsiRK(prev_v_f,v_nzWIdx,v_nzWbinSize,m_RK);	
% 						if strcmp(obj.ch_dataType,'real')
% 							m_IfPsiRK = m_IfPsiRK([1:580,601:2400],:);
% 						end		
% 						v_alpha = ChannelGainMapEstimator.ridgeRegression(m_IfPsiRK',s_check,obj.mu_w,m_K);
% 						
% 						% [S2] estimate v_f_est
% 						v_RKa = m_RK * v_alpha;
% 						v_tempA = zeros(s_gridNum * s_measurementNum,1);
% 						v_tempA(v_nz_PsiRKaIdx) = v_RKa;
% 						m_A = reshape(v_tempA,s_gridNum,s_measurementNum)'; % see the definition of "A" in the paper
% 						if strcmp(obj.ch_dataType,'real')
% 							m_A = m_A([1:580,601:2400],:);
% 						end
% 						
% 						v_f_est = obj.chooseSolver(s_check,m_A');
% 						
% 						estimateDifference = norm(prev_v_f - v_f_est,2);
% 						prev_v_f = v_f_est;
% 						s_iterationNum = s_iterationNum + 1;
% 					end
				case 'simultaneous'
					% we need an estimator jointly estimating v_f_est,
					% pathloss exponent, and tx/rx gains only with shadow-faded channel
					% gain measurements
					s_sensorNum = size(m_sensorPos,2);
					v_sensorDistancesdB = zeros(s_measurementNum,1);
					
					if strcmp(obj.ch_gainType,'different') == 1
						m_Omega = obj.sensorMapOp(m_sensorInd,s_sensorNum);
					else
						m_Omega = ones(s_measurementNum,1);
					end
					
					for s_measurementInd = 1: s_measurementNum
						v_txPos = m_sensorPos(:,m_sensorInd(1,s_measurementInd));
						v_rxPos = m_sensorPos(:,m_sensorInd(2,s_measurementInd));
						v_sensorDistancesdB(s_measurementInd) = 10 * log10(norm(v_txPos-v_rxPos));
					end
					
					m_barI = eye(s_measurementNum) - (1/(norm(v_sensorDistancesdB,2)^2))*(v_sensorDistancesdB)*(v_sensorDistancesdB');
					m_tempOmega = m_barI * m_Omega;
					% add "1e-4 * eye(s_sensorNum)" to "m_tempOmega' * m_tempOmega" for matrix inversion
					m_barOmega = m_Omega * ((m_tempOmega' * m_tempOmega + 1e-4 * eye(size(m_Omega,2))) \ (m_tempOmega' * m_barI));
					m_tildeI = m_barI * (eye(s_measurementNum) - m_barOmega );
					v_barmeasurements = m_tildeI * v_measurements;
					
					s_iterationNum = 1;
					while (estimateDifference > s_stoppingCriterion) && (s_iterationNum < s_iterationMax)
						% [S1] estimate v_alpha
						m_IfR = obj.computeIfR(prev_v_f,v_centroidsInd);
						m_IfRK = m_IfR * m_K;
						m_barIfKR = -1 * m_IfRK' * m_tildeI';
						v_alpha = ChannelGainMapEstimator.ridgeRegression(m_barIfKR,v_barmeasurements,obj.mu_w,m_K);
						
						% [S2] estimate v_f_est
						v_RKa = m_RK * v_alpha;
						m_A = reshape(v_RKa,s_gridNum,s_measurementNum)'; % see the definition of "A" in the paper
						m_barA = -1 * m_A' * m_tildeI';
						v_f_est = obj.chooseSolver(v_barmeasurements,m_barA);
						
						estimateDifference = norm(prev_v_f - v_f_est,2);
						prev_v_f = v_f_est;
						s_iterationNum = s_iterationNum + 1;
					end
					
			end
			
			% blind estimator outputs
			m_F_est = reshape(v_f_est,N_x,N_y);
			h_w_est = @(v_input1,v_input2) obj.represeThm(v_alpha,m_centroids,[v_input1;v_input2]);
			
		end
		
		function [m_centroids, v_centroidsInd,v_alpha,h_w_est] = optfunctionEstimation(obj,t_sensorPos,t_sensorInd)
		% This is a function to estimate optimal kernel coefficients to express a weight function lying in RKHS
		% (accordingly, optimally estimated function, too) using a representer theorem.
		
		% 1. Compute features and correponding kernel matrix. For a
		% practical purpose, features are clustered and represented by
		% much less number of cluster centroids. Those
		% centroids are found by either 'random' sampling, or
		% well-known 'k-means'.
		
		[~, N_y] = size(obj.ini_F);		
		switch obj.ch_dataType
			case 'synthetic'
				m_sensorPos=t_sensorPos;
				m_sensorInd=t_sensorInd;
			case 'real'
				m_sensorPos=t_sensorPos(:,:,2);
				m_sensorInd=t_sensorInd(:,:,2);
		end
			
		[m_centroids, v_centroidsInd]= obj.findCentroids(m_sensorPos,m_sensorInd);
		m_K = ChannelGainMapEstimator.kernelMatrix(m_centroids,obj.h_kernel);


		v_wMeasurements = zeros(size(m_centroids,2),1);
		for s_centroidIdx = 1 : size(m_centroids,2)
			v_wMeasurements(s_centroidIdx,1) = obj.h_w(m_centroids(1,s_centroidIdx),m_centroids(2,s_centroidIdx));
		end
		
		% 2. Estimate v_alpha by solving 
		% 1/(s_centroidsNum) * || v_wMeasurements - m_K * v_alpha ||_2^{2}
		% + \mu_f * v_alpha' * K *v_alpha
		% which is a ridge regression problem. v_wMeasurements are
		% collection of vectorized weight matrices evaluated at all
		% m_centroids).
		
		v_alpha = ChannelGainMapEstimator.ridgeRegression(m_K,v_wMeasurements,obj.mu_w,m_K);
		save('centroids_estimated.mat','m_centroids','v_centroidsInd','v_alpha');
			
		h_w_est = @(v_input1,v_input2) obj.represeThm(v_alpha,m_centroids,[v_input1;v_input2]);	
		
		end

	end
	
	methods
		
		
		
		function v_f_est = chooseSolver(obj,s_check,m_vecWCollection)
			% choose a solver among 'tikhonov'/'l1_ISTA'/'l1_PCO'/'totalvariation'
			% according to a regularization type (ch_reg_f_type).
			[N_x, N_y] = size(obj.ini_F);
			switch obj.ch_reg_f_type
				
				case 'tikhonov'
					switch obj.ch_dataType
						case 'synthetic'
							m_tikhonov = eye(N_x*N_y);
						case 'real'
							m_spatialCov = ChannelGainMapEstimator.spatialCovMat(N_x,N_y,obj.s_resolution);
							m_tikhonov = inv(m_spatialCov);
					end
					v_f_est = ChannelGainMapEstimator.ridgeRegression(m_vecWCollection,s_check,obj.mu_f,m_tikhonov);
				case 'l1_ISTA'
					v_f_est = ChannelGainMapEstimator.lassoISTA(m_vecWCollection,s_check,obj.mu_f,obj.ini_F(:),obj.rho);
				case 'l1_PCO'
					v_f_est = ChannelGainMapEstimator.lassoPCO(m_vecWCollection,s_check,obj.mu_f,obj.ini_F(:));
				case 'totalvariation'
					m_P = ChannelGainMapEstimator.permOpColMaj2RowMaj(N_x,N_y);
					m_D = ChannelGainMapEstimator.differenceOp(obj.ini_F);
					v_f_est = ChannelGainMapEstimator.tv_ADMM(m_vecWCollection,s_check,obj.mu_f,obj.ini_F,obj.rho,m_P,m_D);
			end
			
		end
		
% 		function m_featuresPhi = findPhi(obj,m_sensorPos,m_sensorInd,s_yAxiSize,s_xAxiSize)
% 			% This is a function to find every features only within an ellipse for non-zero weights.
% 			
% 			% OUTPUT:
% 			%   m_featuresPhi    2-by-(s_measurementsNum * s_gridNumInEllipse)
% 			s_measurementNum = size(m_sensorInd,2);
% 			x_axis = repmat((0:s_xAxiSize-1)./obj.s_resolution, s_yAxiSize, 1);   % x_axis of a grid
% 			y_axis = repmat((s_yAxiSize-1:-1:0)'./obj.s_resolution, 1, s_xAxiSize); % y_axis of a grid
% 			s_gridNum = size(x_axis,2) * size(y_axis,2);
% 			cnt = 1;
% 			
% 			m_featuresPhi = zeros(2,s_measurementNum * s_gridNum);
% 			ccnt = 0;
% 			for s_measurementInd = 1 : s_measurementNum
% 				v_txPos = m_sensorPos(:,m_sensorInd(1,s_measurementInd));
% 				v_rxPos = m_sensorPos(:,m_sensorInd(2,s_measurementInd));
% 				
% 				s_phi1 = norm(v_txPos-v_rxPos);
% 				m_phi2 = sqrt( (x_axis-v_txPos(1)).^2 + (y_axis-v_txPos(2)).^2 ) + sqrt( (x_axis-v_rxPos(1)).^2 + (y_axis-v_rxPos(2)).^2 );
% 				v_phi2 = m_phi2(:);
% 				%%%%%%%%%%%%% Check! Do we need to keep all phi's even
% 				%%%%%%%%%%%%% outside of the ellipse?
% 				% 				for s_gridInd = 1 : s_gridNum
% 				% 					if v_phi2(s_gridInd) <= s_phi1 + obj.lambda_W/2
% 				% 						m_featuresPhi(:,cnt) = [s_phi1; v_phi2(s_gridInd)];
% 				% 						cnt = cnt + 1;
% 				% 					end
% 				% 				end
% 				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 				for s_gridInd = 1 : s_gridNum
% 					if v_phi2(s_gridInd) < s_phi1 + obj.lambda_W/2
% 						m_featuresPhi(:,cnt) = [s_phi1; v_phi2(s_gridInd)];
% 						ccnt = ccnt + 1;
% 					else
% 						m_featuresPhi(:,cnt) = [s_phi1; s_phi1 + obj.lambda_W/2];
% 					end
% 					cnt = cnt + 1;
% 				end
% 				
% 			end
% 		end
% 
% 		function [m_centroids, v_centroidsInd]= findCentroids(obj,m_sensorPos,m_sensorInd)
% 			
% 			% OUTPUT:
% 			%   m_centroids       (2-by-s_clusterNum) matrix of centroids
% 			%   m_featuresPhi     (2-by-n_measurements * Ng) matrix of all
% 			%                     feature vectors
% 			%   v_centroidsInd    (n_measurements * Ng)-by-1 vector
% 			%                     of cluster indicies where each feature
% 			%                     vector belongs to
% 			
% 			% Check
% 			clustering = ~isempty(obj.s_clusterNum);
% 			[s_yAxiSize, s_xAxiSize] = size(obj.ini_F);
% 			s_measurementNum = size(m_sensorInd,2);
% 			s_gridNum = s_yAxiSize * s_xAxiSize; % Number of grid points
% 			
% 			
% 			%1. Collect features (phi's) only within the ellipses. A
% 			%function should consider the resolution of the SLF, as well.
% 			m_featuresPhi = obj.findPhi(m_sensorPos,m_sensorInd,s_yAxiSize,s_xAxiSize);
% 			
% 			%2. Find centroids of phi's with the size of s_clusterNum
% 			%according to ch_clustType
% 			
% 			if clustering
% 				
% 				% check
% 				assert( obj.s_clusterNum <= s_gridNum*s_measurementNum );
% 				
% 				switch obj.ch_clustType
% 					case 'kmeans'
% 						% v_centroidsInd := size(m_featuresPhi,2)-by-1 vector
% 						% containing indices of clusters where each column
% 						% of m_collectPhi belongs to
% 						% m_centroids := 2-by-s_clusterNum matrix where each
% 						% column is a centroid of a cluster
% 						[v_centroidsInd,m_centroids] = kmeans(m_featuresPhi',obj.s_clusterNum,'Replicates',2,'MaxIter',500,'start','sample','emptyaction','singleton');
% 						
% 						if  (sum(sum(isnan(m_centroids)))>0) || ( sum(sum(m_centroids==Inf))>0 )
% 							error('kmeans returned NaN');
% 						end
% 						[m_centroids, ia, ic] = unique(m_centroids,'rows'); % Rmove redundant centroids
% 						v_centroidsInd = ic(v_centroidsInd);
% 						obj.s_clusterNum = size(ia,1);
% 						m_centroids = m_centroids';
% 					case 'random'
% 						% choose obj.s_clusterNum feature vectors (cols of phi_col) --> centroids
% 						m_tempCollectPhi = m_featuresPhi';
% 						m_UniqcollectPhi = unique(m_tempCollectPhi,'rows')';
% 						
% 						v_rndCentroidInd = randperm(size(m_UniqcollectPhi,2),obj.s_clusterNum);
% 						m_centroids = m_UniqcollectPhi(:,v_rndCentroidInd);
% 						
% 						m_Dist = pdist2(m_featuresPhi',m_centroids');% matD := distance matrix between coordinates and centroids
% 						[~,v_centroidsInd] = min(m_Dist,[],2);
% 				end
% 				
% 				
% 			else  % no clustering
% 				m_centroids = m_featuresPhi; % evl_pnt := evluation point of kerenl.
% 				v_centroidsInd = 1:1:size(m_featuresPhi,2);
% 			end
% 			
% 			% 			ChannelGainMapEstimator.plot_clusters(m_centroids,m_featuresPhi,v_centroidsInd);
% 			
% 		end
		
		function [m_featuresPhi,m_featuresPhi4Centroids] = findPhi(obj,m_sensorPos,m_sensorInd,s_yAxiSize,s_xAxiSize)
			% This is a function to find every features only within an ellipse for non-zero weights.
			
			% OUTPUT:
			%   m_featuresPhi    2-by-(s_measurementsNum * s_gridNumInEllipse)
			s_measurementNum = size(m_sensorInd,2);
			x_axis = repmat((0:s_xAxiSize-1)./obj.s_resolution, s_yAxiSize, 1);   % x_axis of a grid
			y_axis = repmat((s_yAxiSize-1:-1:0)'./obj.s_resolution, 1, s_xAxiSize); % y_axis of a grid
			s_gridNum = size(x_axis,2) * size(y_axis,2);
			cnt = 1;
			
			m_featuresPhi = zeros(2,s_measurementNum * s_gridNum);
			ccnt = 1;
			for s_measurementInd = 1 : s_measurementNum
				v_txPos = m_sensorPos(:,m_sensorInd(1,s_measurementInd));
				v_rxPos = m_sensorPos(:,m_sensorInd(2,s_measurementInd));
				
				s_phi1 = norm(v_txPos-v_rxPos);
				m_phi2 = sqrt( (x_axis-v_txPos(1)).^2 + (y_axis-v_txPos(2)).^2 ) + sqrt( (x_axis-v_rxPos(1)).^2 + (y_axis-v_rxPos(2)).^2 );
				v_phi2 = m_phi2(:);
				
				for s_gridInd = 1 : s_gridNum
					m_featuresPhi(:,cnt) = [s_phi1; v_phi2(s_gridInd)];
					
					if v_phi2(s_gridInd) < s_phi1 + obj.lambda_W/2
						m_featuresPhi4Centroids(:,ccnt) = [s_phi1; v_phi2(s_gridInd)];
						ccnt = ccnt + 1;
					end
					cnt = cnt + 1;						
				end
				
				m_featuresPhi4Centroids(:,ccnt) = [s_phi1; s_phi1 + obj.lambda_W/2];
				ccnt = ccnt + 1;
			end
		end
		
		function [m_centroids, v_centroidsInd]= findCentroids(obj,m_sensorPos,m_sensorInd)
			
			% OUTPUT:
			%   m_centroids       (2-by-s_clusterNum) matrix of centroids
			%   m_featuresPhi     (2-by-n_measurements * Ng) matrix of all
			%                     feature vectors
			%   v_centroidsInd    (n_measurements * Ng)-by-1 vector
			%                     of cluster indicies where each feature
			%                     vector belongs to
			
			% Check
			clustering = ~isempty(obj.s_clusterNum);
			[s_yAxiSize, s_xAxiSize] = size(obj.ini_F);
			s_measurementNum = size(m_sensorInd,2);
			s_gridNum = s_yAxiSize * s_xAxiSize; % Number of grid points
			
			
			%1. Collect features (phi's) only within the ellipses. A
			%function should consider the resolution of the SLF, as well.
			[m_featuresPhi,m_featuresPhi4Centroids] = obj.findPhi(m_sensorPos,m_sensorInd,s_yAxiSize,s_xAxiSize);
			
			%2. Find centroids of phi's with the size of s_clusterNum
			%according to ch_clustType
			
			if clustering
				
				% check
				assert( obj.s_clusterNum <= s_gridNum*s_measurementNum );
				
				switch obj.ch_clustType
					case 'kmeans'
						% v_centroidsInd := size(m_featuresPhi,2)-by-1 vector
						% containing indices of clusters where each column
						% of m_collectPhi belongs to
						% m_centroids := 2-by-s_clusterNum matrix where each
						% column is a centroid of a cluster
						[~,m_centroids] = kmeans(m_featuresPhi4Centroids',obj.s_clusterNum,'Replicates',2,'MaxIter',500,'start','sample','emptyaction','singleton');
						
						if  (sum(sum(isnan(m_centroids)))>0) || ( sum(sum(m_centroids==Inf))>0 )
							error('kmeans returned NaN');
						end
						[m_centroids, ia, ic] = unique(m_centroids,'rows'); % Rmove redundant centroids
						m_Dist = pdist2(m_featuresPhi',m_centroids);% matD := distance matrix between coordinates and centroids
						[~,v_centroidsInd] = min(m_Dist,[],2);
						m_centroids = m_centroids';
					case 'random'
						% choose obj.s_clusterNum feature vectors (cols of phi_col) --> centroids
						m_tempCollectPhi = m_featuresPhi4Centroids';
						m_UniqcollectPhi = unique(m_tempCollectPhi,'rows')';
						
						v_rndCentroidInd = randperm(size(m_UniqcollectPhi,2),obj.s_clusterNum);
						m_centroids = m_UniqcollectPhi(:,v_rndCentroidInd);
						m_Dist = pdist2(m_featuresPhi',m_centroids');% matD := distance matrix between coordinates and centroids
						[~,v_centroidsInd] = min(m_Dist,[],2);
				end
				
				
			else  % no clustering
				m_centroids = m_featuresPhi; % evl_pnt := evluation point of kerenl.
				v_centroidsInd = 1:1:size(m_featuresPhi,2);
			end
			
			ChannelGainMapEstimator.plot_clusters(m_centroids,m_featuresPhi,v_centroidsInd);
			
		end
			
		function [m_featuresPhi,m_featuresPhi4Centroids,v_nzWIdx,v_nzWbinSize] = findPhiEllip(obj,m_sensorPos,m_sensorInd,s_yAxiSize,s_xAxiSize)
			% This is a function to find features only within an ellipse for non-zero weights.
			
			% OUTPUT:
			%   m_featuresPhi    2-by-(s_measurementsNum * s_gridNumInEllipse)
			%   v_nzWIdx         (sum of L_t)-by-1 vector contains indices
			%                    corresponding to nonzero weight, where L_t is the number of
			%                    nonzero weight at time t.
			%   v_nzWbinSize     (s_measurementsNum-by-1) each element is the number of
			%                    nonzero weights for each time slot.
			
			s_measurementNum = size(m_sensorInd,2);
			x_axis = repmat((0:s_xAxiSize-1)./obj.s_resolution, s_yAxiSize, 1);   % x_axis of a grid
			y_axis = repmat((s_yAxiSize-1:-1:0)'./obj.s_resolution, 1, s_xAxiSize); % y_axis of a grid
			s_gridNum = size(x_axis,2) * size(y_axis,2);
			s_featuresNum = 1;
			s_featuresCentNum = 1;
			for s_measurementInd = 1 : s_measurementNum
				v_txPos = m_sensorPos(:,m_sensorInd(1,s_measurementInd));
				v_rxPos = m_sensorPos(:,m_sensorInd(2,s_measurementInd));
				
				s_phi1 = norm(v_txPos-v_rxPos);
				m_phi2 = sqrt( (x_axis-v_txPos(1)).^2 + (y_axis-v_txPos(2)).^2 ) + sqrt( (x_axis-v_rxPos(1)).^2 + (y_axis-v_rxPos(2)).^2 );
				v_phi2 = m_phi2(:);
				s_nzWbinSize = 0;
				for s_gridInd = 1 : s_gridNum
					
					if v_phi2(s_gridInd) < s_phi1 + obj.s_SemiAxisLength4Sample /2
						m_featuresPhi(:,s_featuresNum) = [s_phi1; v_phi2(s_gridInd)];
						v_nzWIdx(s_featuresNum,1) = s_gridInd;
						s_featuresNum = s_featuresNum + 1;
						s_nzWbinSize = s_nzWbinSize + 1;
						
					end
					
					if v_phi2(s_gridInd) < s_phi1 + obj.lambda_W/2
						m_featuresPhi4Centroids(:,s_featuresCentNum) = [s_phi1; v_phi2(s_gridInd)];
						s_featuresCentNum = s_featuresCentNum + 1;
					end
					
					
				end
				v_nzWbinSize(s_measurementInd,1) = s_nzWbinSize;
				
			end
		end

		function [m_centroids, v_centroidsInd, v_nzWIdx,v_nzWbinSize]= findCentroidsEllip(obj,m_sensorPos,m_sensorInd)
			
			% OUTPUT:
			%   m_centroids       (2-by-s_clusterNum) matrix of centroids
			%   m_featuresPhi     (2-by-n_measurements * Ng) matrix of all
			%                     feature vectors
			%   v_centroidsInd    (n_measurements * Ng)-by-1 vector
			%                     of cluster indicies where each feature
			%                     vector belongs to
			%   v_nzWIdx         (sum of L_t)-by-1 vector contains indices
			%                    corresponding to nonzero weight, where L_t is the number of
			%                    nonzero weights at time t.
			%   v_nzWbinSize     (s_measurementsNum-by-1) each element is the number of
			%                    nonzero weights for each time slot.
			
			% Check
			clustering = ~isempty(obj.s_clusterNum);
			[s_yAxiSize, s_xAxiSize] = size(obj.ini_F);
			s_measurementNum = size(m_sensorInd,2);
			s_gridNum = s_yAxiSize * s_xAxiSize; % Number of grid points
			
			
			%1. Collect features (phi's) only within the ellipses. A
			%function should consider the resolution of the SLF, as well.
			[m_featuresPhi,m_featuresPhi4Centroids,v_nzWIdx,v_nzWbinSize] = obj.findPhiEllip(m_sensorPos,m_sensorInd,s_yAxiSize,s_xAxiSize);
			
			%2. Find centroids of phi's with the size of s_clusterNum
			%according to ch_clustType
			
			if clustering
				
				% check
				assert( obj.s_clusterNum <= s_gridNum*s_measurementNum );
				
				switch obj.ch_clustType
					case 'kmeans'
						% v_centroidsInd := size(m_featuresPhi,2)-by-1 vector
						% containing indices of clusters where each column
						% of m_collectPhi belongs to
						% m_centroids := 2-by-s_clusterNum matrix where each
						% column is a centroid of a cluster
						[~,m_centroids] = kmeans(m_featuresPhi4Centroids',obj.s_clusterNum,'Replicates',2,'MaxIter',500,'start','sample','emptyaction','singleton');
						
						if  (sum(sum(isnan(m_centroids)))>0) || ( sum(sum(m_centroids==Inf))>0 )
							error('kmeans returned NaN');
						end
						[m_centroids, ia, ic] = unique(m_centroids,'rows'); % Rmove redundant centroids
						m_Dist = pdist2(m_featuresPhi',m_centroids);% matD := distance matrix between coordinates and centroids
						[~,v_centroidsInd] = min(m_Dist,[],2);
						m_centroids = m_centroids';
					case 'random'
						% choose obj.s_clusterNum feature vectors (cols of phi_col) --> centroids
						m_tempCollectPhi = m_featuresPhi4Centroids';
						m_UniqcollectPhi = unique(m_tempCollectPhi,'rows')';
						
						v_rndCentroidInd = randperm(size(m_UniqcollectPhi,2),obj.s_clusterNum);
						m_centroids = m_UniqcollectPhi(:,v_rndCentroidInd);
						
						m_Dist = pdist2(m_featuresPhi',m_centroids');% matD := distance matrix between coordinates and centroids
						[~,v_centroidsInd] = min(m_Dist,[],2);
				end
				
				
			else  % no clustering
				m_centroids = m_featuresPhi; % evl_pnt := evluation point of kerenl.
				v_centroidsInd = [1:1:size(m_featuresPhi,2)]';
			end
			
			% 			ChannelGainMapEstimator.plot_clusters(m_centroids,m_featuresPhi,v_centroidsInd);
			
		end
		
		function m_IfR = computeIfR(obj,v_f,v_centroidsInd)
			% This is a function to efficiently calculate "kron(eye(s_measurementNum),v_f')*m_R".
			
			% fR_mat = ( I \kron f^T ) R
			% i) R_tilde:=[R_1' ; ... ; R_T'] * f
			% ii) Rf_vec := sum(R_tilde,2)
			% iii) fR_Mat = unvec(Rf_vec)'
			
			s_numGrid = size(v_f,1);
			s_measurementNum = length(v_centroidsInd)/ s_numGrid;
			idx_phi_mat = reshape(v_centroidsInd,s_numGrid,s_measurementNum); % non-zero locations of the columns of R matrix
			row_inds = idx_phi_mat + obj.s_clusterNum*ones(s_numGrid,1)*(0:s_measurementNum-1);
			col_inds = repmat(1:s_numGrid,1,s_measurementNum)';
			entry_vals = repmat(v_f,s_measurementNum,1);
			
			v_Rf = full(sum(sparse(row_inds(:),col_inds,entry_vals,obj.s_clusterNum*s_measurementNum,s_numGrid,length(v_centroidsInd)),2));
			
			m_IfR = reshape(v_Rf,obj.s_clusterNum,s_measurementNum)';
			
		end
	
		function m_IfPsiRK = computeIfPsiRK(obj,v_f,v_nzWIdx,v_nzWbinSize,m_RK)
			% This is a function to efficiently calculate "kron(eye(s_measurementNum),v_f'*Phi_i)*m_R * m_K".
			
			m_pf = repmat(v_f(v_nzWIdx),1,obj.s_clusterNum);
			m_tempMat = m_pf .* m_RK;
			m_IfPsiRK(1,:) = sum(m_tempMat(1:v_nzWbinSize(1,1),:));
			for s_ItrIndex = 2 : size(v_nzWbinSize,1)
				s_lastPoint = sum(v_nzWbinSize(1:(s_ItrIndex-1),1));
				m_IfPsiRK(s_ItrIndex,:) = sum(m_tempMat(s_lastPoint+1:s_lastPoint + v_nzWbinSize(s_ItrIndex,1),:));
			end

		end

		
		function h_func = represeThm(obj,v_alpha,m_featureVecs,v_input)
			% Learn a function by using the representer theorem.
			
			% INPUT:
			%   v_alpha             (t*Ng)-by-1 coefficients of kernel
			%   m_featureVecs       2-by-(t*Ng) matrix of feature vectors
			%                       corresponding to each coefficient
			%   v_input             2-by-1 kernel input
			%   h_kernel            Kernel function
			% OUTPUT:
			%   h_func              estimated function
			
			s_coefficientNum = size(v_alpha,1);
			
			v_K = zeros(1,s_coefficientNum);
			for i = 1 : s_coefficientNum
				v_K(1,i) = obj.h_kernel(v_input,m_featureVecs(:,i));
			end
			% 			v_concateinput = kron(ones(1,s_coefficientNum),v_input);
			% 			v_K = obj.h_kernel(v_concateinput,m_featureVecs);
			
			h_func = v_K * v_alpha;
		end
		
	end


    methods(Static) % solvers

        % COMMON INPUT:
        %      m_A         (N_g-by-s_measurementNum) parameter matrix 
        %      v_s         s_measurementNum-by-1 vector with (noisy/noiseless) measurements 
        %      s_mu_f      a regularization weight

        % COMMON OUTPUT:
        %      v_f         (N_g-by-1) (vectorized) spatial loss field
        
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1. Ridge Regression 
        % min_{v_f} (1/s_measurementNum) * || v_s - m_A' * v_f ||_2^{2} +
        % s_mu_f * v_f' * m_tikhonov * v_f
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function v_f = ridgeRegression(m_A,v_s,s_mu_f,m_tikhonov)
            [s_sizex,s_measurementNum] = size(m_A);
            v_f = (m_A * m_A' + s_measurementNum * m_tikhonov * s_mu_f)\(m_A *v_s);
        end
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2. LASSO via Iterative Shrinkage/Thresholding Algorithm (ISTA)
        % min_{v_f} (1/s_measurementNum) * || v_s - m_A' * v_f ||_2^{2} + s_mu_f * ||v_f||_1
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function v_f = lassoISTA(m_A,v_s,s_mu_f,ini_v_f,rho)
            % INPUT:
            %    ini_v_f          (N_g-by-1) initialization of a spatial loss field
            %    rho              a step size for a gradient descent 

            s_measurementNum = size(m_A,2);
            
            % Simulation setup            
            prev_v_f = ini_v_f;
            s_stoppingCriterion = 1e-5;
            estimateDifference = inf;
            
            s_iterationNum = 1;
            while (estimateDifference > s_stoppingCriterion) && (s_iterationNum < 1e4)

                % [Step1] Gradient descent
                grad = -2 * m_A *(v_s - m_A' * prev_v_f)/s_measurementNum;
                v_z = prev_v_f - rho * grad;
                
                % [Step2] Proximal operation
                v_f = wthresh(v_z,'s', rho*s_mu_f);
                
                estimateDifference = norm(prev_v_f - v_f,2);
                prev_v_f = v_f;
                s_iterationNum = s_iterationNum + 1;
            end
            

        end
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 3. LASSO via Pathwise Coordinate Optimization (PCO)
        % min_{v_f} (1/s_measurementNum)*|| v_s - m_A' * v_f ||_2^{2} + s_mu_f * ||v_f||_1
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function v_f = lassoPCO(m_A,v_s,s_mu_f,ini_v_f)
            % INPUT:
            %    ini_v_f          (N_g-by-1) initialization of a spatial loss field

            [N_g,s_measurementNum] = size(m_A);
            
            % Simulation setup            
            prev_v_f = ini_v_f;
            s_stoppingCriterion = 1e-5;
            estimateDifference = inf;
                                 
            s_iterationNum = 1;
            while (estimateDifference > s_stoppingCriterion) && (s_iterationNum < 1e4)
                
                % Precomputation for sped up
                v_temp_tilde_s = v_s - m_A' * prev_v_f;
                
                for l = 1 : N_g
                    v_a_l = m_A(l,:)';
                    v_tilde_s = v_temp_tilde_s + v_a_l * prev_v_f(l,1);
                    if norm(v_a_l,2) ~= 0
                        v_f(l,1) = wthresh(v_tilde_s' * v_a_l,'s', s_measurementNum * s_mu_f * 0.5)/(norm(v_a_l,2)^2);
                    else
                        v_f(l,1) = 0;
                    end
                    v_temp_tilde_s = v_tilde_s - v_a_l * v_f(l,1);
                end
                              
                estimateDifference = norm(prev_v_f - v_f,2);
                prev_v_f = v_f;
                s_iterationNum = s_iterationNum + 1;
            end

        end
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 4. 2D-Total Variation via Alternating direction method of multipliers(ADMM)
        % min_{v_f} (1/s_measurementNum)*|| v_s - m_A' * v_f ||_2^{2} + s_mu_f * ||v_f||_TV
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function v_f = tv_ADMM(m_A,v_s,s_mu_f,ini_m_F,rho,m_P,m_D)
            % INPUT:
            %    ini_v_f          (N_g-by-1) initialization of a spatial loss field
            %    rho              a step size for a gradient descent 
            %    m_P              (N_g-by-N_g) permutation operator  
            %    m_D              (N_y*(N_x-1) + N_x*(N_y-1))-by-N_g
            %                     2-dimensional difference operator
                      
            s_measurementNum = size(m_A,2);
            [N_x,N_y] = size(ini_m_F);
            
            % Simulation setup
            prev_v_f = ini_m_F(:);
            s_stoppingCriterion = 1e-4;
            estimateDifference = inf;
            
            m_D1 = m_D(1:N_y*(N_x-1),:) ;
            m_D2 = m_D(N_y*(N_x-1) + 1 : N_y*(N_x-1) + N_x*(N_y-1),:);
            % precomputation for sped up
            til_A = (2/s_measurementNum) * m_A * m_A';
            til_D = m_D1'*m_D1;
            DP = m_D2 * m_P;
            til_P = DP' * DP;
            til_H = (til_A + rho * (til_D + til_P) );
            til_s = (2/s_measurementNum) * m_A * v_s;
            
            % Initialization of dual/auxiliary variables
            v_gamma_x = zeros(size(m_D1,1),1); % dual variable
            v_gamma_y = zeros(size(m_D2,1),1); % dual variable
            v_d_x = zeros(size(m_D1,1),1); % aux. variable v_d_x = m_D1 * v_f
            v_d_y = zeros(size(m_D2,1),1); % aux. variable v_d_y = m_D2 * m_P * v_f
                        
            s_iterationNum = 1;
            while (estimateDifference > s_stoppingCriterion) && (s_iterationNum < 1 * 1e2)

                % [S1] Dual ascent for dual variables
                v_gamma_x = v_gamma_x + rho * (m_D1*prev_v_f - v_d_x);
                v_gamma_y = v_gamma_y + rho * (DP*prev_v_f - v_d_y);
                % [S2] Updates of aux variables d_x and d_y
                v_d_x = wthresh(m_D1*prev_v_f + v_gamma_x./rho,'s', s_mu_f/rho);
                v_d_y = wthresh(DP*prev_v_f + v_gamma_y./rho,'s', s_mu_f/rho);

                % [S3] Updates of SLF
                v_f = ( til_H ) \(m_D1' *(rho * v_d_x - v_gamma_x) + DP' *(rho * v_d_y -v_gamma_y) + til_s); 
                                            
                estimateDifference = norm(prev_v_f - v_f,2);
                prev_v_f = v_f;
                s_iterationNum = s_iterationNum + 1;
            end
        end
        
		
    end
    
    methods(Static) % Utilities
        function m_P = permOpColMaj2RowMaj(N_x,N_y)
        % Generate a permuation operator to change a column-major form vectorization of a matrix to
        % row-major form, i.e., m_P * m_x = m_X'(:) where m_x = m_X(:)
        
        % Output:
        %   m_P            N_g-by-N_g permutation operator  

            N_g = N_x * N_y;
            v_seed = [1,zeros(1,N_g-1)];
            for i = 1 : N_y
                m_seed(i,:) = circshift(v_seed,[0,(i-1)*N_x]);
            end

            m_P = m_seed;
            for i = 2 : N_x
                m_P = [m_P;circshift(m_seed,[0,(i-1)])];
            end

        end
        
        function m_D = differenceOp(m_X)
        % Generate a 2D difference operator to give ||m_X||_{2D-TV anistropic} = ||m_D * v_x||_1
        % Output:
        %    m_D        (N_y*(N_x-1) + N_x*(N_y-1))-by-N_g 0,+-1 matrix where  
        %                m_D = [m_D1; m_D2], dim(m_D1) = N_y*(N_x-1)-by-Ng, dim(m_D2) = N_x*(N_y-1)-by-Ng               
      
            [N_x,N_y] = size(m_X);
            N_g = N_x * N_y;

            % 1-D
            v_seed1 = [1,-1,zeros(1,N_x-2)];
            m_seed1 = zeros(N_x-1,N_x);
            for i = 1 : N_x-1
                m_seed1(i,:) = circshift(v_seed1,[1,i-1]);
            end
            m_D1 = kron(eye(N_y),m_seed1);

            % 2-D
			v_seed2 = [1,-1,zeros(1,N_y-2)];
			m_seed2 = zeros(N_y-1,N_y);
			for i = 1 : N_y-1
				m_seed2(i,:) = circshift(v_seed2,[1,i-1]);
			end
			m_D2 = kron(eye(N_x),m_seed2);
		
            m_D = [m_D1;m_D2];
		end
  
		function m_Omega = sensorMapOp(m_sensorInd,s_sensorNum)
        % Generate a (s_measurementNum-by-s_sensorNum) matrix to map a pair of sensors on each measurement.
		% Each row contains two 1s on the locations of selected sensors, and 0s on the others.
                   
            s_measurementNum = size(m_sensorInd,2);
			m_Omega = zeros(s_measurementNum,s_sensorNum);
			for s_measurementIdx = 1 : s_measurementNum
				m_Omega(s_measurementIdx,m_sensorInd(:,s_measurementIdx)')=1;
			end

		end
		
		function m_spatialCov = spatialCovMat(s_yAxiSize,s_xAxiSize,s_resolution)
			% Generate a (s_xAxiSize*s_yAxiSize-by-s_xAxiSize*s_yAxiSize) spatial covariance
			% matrix. Covariance between two points v_x and v_y is "exp(-1
			% * norm(v_x - v_y) )".
				
			x_axis = repmat((0:s_xAxiSize-1)./s_resolution, s_yAxiSize, 1);   % x_axis of a grid
			y_axis = repmat((s_yAxiSize-1:-1:0)'./s_resolution, 1, s_xAxiSize); % y_axis of a grid
					
			v_gridY = y_axis(:);
			v_girdX = x_axis(:);
			
			s_gridYcovMat = size(y_axis,2);
			s_gridXcovMat = size(x_axis,2);
			
			u = 1;
			for i = 1: s_gridYcovMat * s_gridXcovMat
				for j = 1 : s_gridYcovMat * s_gridXcovMat
					v_spatialCov(u,1) = sqrt((v_girdX(j,1) - v_girdX(i,1))^2 + (v_gridY(j,1) - v_gridY(i,1))^2);
					u = u + 1;
				end
			end
			m_temp_spatialCov  = reshape(v_spatialCov,s_gridYcovMat*s_gridXcovMat,s_gridYcovMat*s_gridXcovMat);
			m_spatialCov = exp(-1*m_temp_spatialCov);
			
		end
				
		function m_K = kernelMatrix(m_centroids,h_kernel)
			% Kernel matrix generator
			
			% INPUT:
			%   m_centroids          2-by-s_centroidNum matrix of centroids
			%   h_kernel             Kernel function
			%
			% OUTPUT:
			%       K               s_centroidNum x s_clusterNum kernel matrix
			
			s_centroidNum = size(m_centroids,2);
			m_K = zeros(s_centroidNum);
			
			for m0 = 1 : s_centroidNum
				for m1 = 1 : s_centroidNum
					m_K(m0,m1) = h_kernel(m_centroids(:,m0),  m_centroids(:,m1) );
				end
			end
			
			min_eig = min(eig(m_K));
			if min_eig <= 0
				m_K = m_K - min_eig* eye(size(m_K,1));
			end
% 			m_K = m_K + 10^8*eps* eye(s_centroidNum);
			
			if rank(m_K) < s_centroidNum
				s_centroidNum
				rank_K = rank(m_K)
			end
					
		end

		function v_nz_PsiRKaIdx = nzRelativeIdxFinder(v_nzWIdx,v_nzWbinSize,s_gridNum)
			% This is a function to find relative nonzero indices for 
			% a aggregated matrix consisted of {Phi_i^{T} * R_{i} * m_kernel
			% * v_alpha}_{i=1}^{s_measurementNum}.
			
			% INPUT:
			%   m_centroids          2-by-s_centroidNum matrix of centroids
			%   h_kernel             Kernel function
			%
			% OUTPUT:
			%   v_nz_PsiRKaIdx         (s_measurementNum * s_gridNum) matrix
			%                        contains relative positions of nonzero
			%                        weights corresponding to R * m_Kernel
			%                        * v_alpha.
			
			v_nz_PsiRKaIdx = v_nzWIdx(1:v_nzWbinSize(1),1);
			for s_ItrIndex = 1 : size(v_nzWbinSize,1)-1
				s_lastPoint = sum(v_nzWbinSize(1:s_ItrIndex,1));
				v_nz_PsiRKaIdx =[ v_nz_PsiRKaIdx; v_nzWIdx(s_lastPoint+1 : s_lastPoint + v_nzWbinSize(s_ItrIndex+1),1) + s_ItrIndex * s_gridNum];		
			end
			
		end

	end
	
	methods(Static) 
		
		function [v_gains_est,s_pathLossExponent_est] = estimateSensorGainAndPathLoss(m_sensorPos,m_sensorInd,v_measurementsNoShadowing,m_Omega)
			% Estimate(or assign) the sumOfGain and pathloss for every pair
			% of sensors when obj.ch_calibrationType is either "none", or
			% "previous".
			
			% m_Omega is optional. See properties of this class.
			
			s_measurementNum = size(m_sensorInd,2);
			%s_sensorNum = size(m_sensorPos,2);
			v_sensorDistances = zeros(s_measurementNum,1);
			
			if (nargin<4) || isempty(m_Omega)
				%strcmp(obj.ch_gainType,'different') == 1
				%m_Omega = obj.sensorMapOp(m_sensorInd,s_sensorNum);
				m_Omega = ones(s_measurementNum,1);
			end
			
			for s_measurementInd = 1: s_measurementNum
				v_txPos = m_sensorPos(:,m_sensorInd(1,s_measurementInd));
				v_rxPos = m_sensorPos(:,m_sensorInd(2,s_measurementInd));
				v_sensorDistances(s_measurementInd) = 10 * log10(norm(v_txPos-v_rxPos));
			end
			
			if  isempty(v_measurementsNoShadowing) 
				error();
			else % 
				% v_estGains / s_estPathLossExponent should be estimated
				m_regMat = [m_Omega,-1 * v_sensorDistances];
				v_parameters = (m_regMat'*m_regMat + 1e-12* eye(size(m_Omega,2) + 1))\(m_regMat'*v_measurementsNoShadowing);
				v_gains_est = v_parameters(1:size(m_Omega,2),1);
				s_pathLossExponent_est = v_parameters(size(m_Omega,2)+1,1);
			end
			
% 			% v_estGains and s_estPathLossExponent should be known or
% 			% estimated before this step.
% 			
% 			v_sumOfGains = m_Omega * v_gains_est;
% 			v_pathLoss = s_pathLossExponent_est * v_sensorDistances;
		end
	
		function m_imageOut = postprocess(m_imageIn)
			
			m_imageOut = m_imageIn;
			m_imageOut( m_imageOut < 0 ) = 0;
			m_imageOut( m_imageOut > 1 ) = 1;
			
		end

		function m_imageOut = postprocessfunction(m_imageIn,m_TrueimageIn,s_mu_out)
			% This is a output image postprocess function to give a
			% visually pleasing output image.
			[N_x,N_y] = size(m_imageIn);
			v_f = m_imageIn(:);
			v_gndTruthf = m_TrueimageIn(:);
			s_gridNum = N_x * N_y;
			cvx_begin
			variable s_a
			variable s_b
			variable v_e(s_gridNum,1)
			minimize ( sum_square(v_gndTruthf - s_a * v_f - s_b * ones(s_gridNum,1) + v_e) + s_mu_out * norm(v_e))
			cvx_end
			
			v_imageOut = s_a * v_f + s_b * ones(s_gridNum,1);
			m_imageOut = reshape(v_imageOut,N_x,N_y);
			
		end
		
		function m_imageOut = postprocessReal(m_imageIn)
			
			m_imageOut = m_imageIn;
			m_imageOut( m_imageOut < 0 ) = 0;
			m_imageOut( m_imageOut > 0.3 ) = 0.3;
			
		end
		
		function plot_clusters(m_centroids,m_featuresPhi,v_centroidsInd)
			% Plot phi's and clusters
			figure;
			plot(m_centroids(1,:),m_centroids(2,:),'rx')
			hold all
			for i = 1 : size(m_centroids,2)
				%phi = phi';
				plot(m_featuresPhi(1,v_centroidsInd==i),m_featuresPhi(2,v_centroidsInd==i),'.')
			end
			
			xlabel('\phi_1')
			ylabel('\phi_2')			
		end
		
		function [m_w_o,m_w_hat] = evaluate_w(h_w,h_w_est,v_rangPhi1,v_rangPhi2,v_intervGrid)
			%  Return functions values of h_w and h_est_w evaluated within
			%  v_Range.
			%
			% INPUT:
			%   w             Original function
			%   v_rangPhi1    (2-by-1) First entry is the minimum value of a range for Phi1, and
			%                 the second entry is the maximum value of a range for Phi1.
			%   v_rangPhi2    (2-by-1) First entry is the minimum value of a range for Phi2, and
			%                 the second entry is the maximum value of a range for Phi2.
			%   v_intervGrid  (2-by-1) interval between grid points over
			%                 Phi1 and Phi2 axis. First entry is of Phi1
			%                 and the second entry is of Phi2.
			
			% OUTPUT:
			%   w_o           Evaluation of ground truth h_w
			%   w_hat         Evaluation of estimated function h_est_w
					
			rng_phi1 = v_rangPhi1(1):v_intervGrid(1):v_rangPhi1(2);
			rng_phi2 = v_rangPhi2(1):v_intervGrid(2):v_rangPhi2(2);
			
			len_phi1 = size(rng_phi1,2);
			len_phi2 = size(rng_phi2,2);
			for i = 1 : len_phi1
				for j = 1 : len_phi2
					if rng_phi1(i) > rng_phi2(j)
						m_w_hat(i,j) = NaN;
						m_w_o(i,j) = NaN;
					else
						m_w_hat(i,j) = h_w_est(rng_phi1(i), rng_phi2(j));
						m_w_o(i,j) = h_w(rng_phi1(i), rng_phi2(j));					
					end
				end
			end
			
					
		end
		
		
	end
	
end

