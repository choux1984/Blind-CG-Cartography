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
		m_tikhonov	 % tikhonov matrix to define a tikhonov regularizer, i.e., R(v_f) = v_f' * m_tik * v_f 
		mu_f         % Regularization parameter for spatial loss field
				
		ch_estimationType
		             % 'non-blind'
					 % 'blind'
					 %
		
        % calibration
		ch_calibrationType = 'none';
		             %  'none'  : values of pathloss and gains taken from
		             %            the properties below
					 %  'previous' : first calibration and then estimation
					 %            of the SLF
					 %  'simultaneous' 
		s_pathLossExponent; 
		v_gains;
					 
		% estimation
		rho          % step size for (ISTA/ADMM) algorithm
		s_resolution = 1;% scaling factor to decide a resolution of an output image, e.g s_resolution = 1 (default)
		                % resolution of the output changes only for real dataset.              
		 
		% non-blind estimation
		h_w = [];    %  
		
		% blind estimation
		mu_w         % Regularization parameter for weight function		
		h_kernel;    % Kernel function to interpolate a weight function. e.g. Gaussian kernel
		s_clusterNum % Number of clusters; Set Nc = [] for no clustering
		ch_clustType % Types of feature clusterting, and corresponding cluster centroids. 
		             % 'random' : random sampling 
					 % 'kmeans' : k-means
		lambda_W     % h_w is estimated with domain equal to an ellipse whose semi-minor axis is lambda_W
	end
		
	methods % estimator
		
		function obj = ChannelGainMapEstimator(varargin)
			obj@Parameter(varargin{:});
		end
				
		function [m_F_est,h_w_est,m_centroids] = estimate(obj,m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing)	
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
					
				case 'blind'
					[m_F_est,h_w_est] = obj.blindEstimation(m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing);
					% (finish this)
					
					
					
					
% 					assert( obj.Nc <= Ng*t );
% 					[evl_pnt,idx_phi,phi_col]  = generate_feature_vecs(N_x,N_y,Tx_pos,Rx_pos,Nc,clustering_type,lambda_W);
% 					K = myKmatrix(evl_pnt,myKfunc);
% 					
% 					[alpha,est_f] = blind_est(K,idx_phi,s_check,mu_w,ini_F(:),eps_error,max_itr,mu_f,reg_f_type,rho);
% 					w = @(input)  myRepThm(alpha,evl_pnt,input,myKfunc);
% 					est_F = reshape(est_f,N_x,N_y);
% 				
									
				otherwise
					error('unrecognized option');
			end
				 
		end
		
		function m_F_est  = nonblindEstimation(obj,m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing)
			% write this using obj.ch_reg_f_type, obj.h_w and obj.mu_f
			% Estimate a spatial loss field in non-blind fashion with a
			% chosen regularizer (tikhonov, l-1, and total variation).
            
			% check
			if isempty(obj.ini_F)&&(strcmp(obj.ch_reg_f_type,'tikhonov')==0)
				assert(~isempty(obj.ini_F));				
            end
            
            % 1. Compute a weight matrix for each pair of Tx / Rx from m_txPos / m_rxPos
            [N_x, N_y] = size(obj.ini_F);
			s_gridNum = N_x * N_y;
            s_measurementNum = length(v_measurements);			
			x_axis = repmat((0:N_y-1)./obj.s_resolution, N_x, 1);   % x_axis of a grid
			y_axis = repmat((N_x-1:-1:0)'./obj.s_resolution, 1, N_y); % y_axis of a grid

            m_vecWCollection = zeros(s_gridNum,s_measurementNum);
            for s_measurementInd = 1 : s_measurementNum
                v_txPos = m_sensorPos(:,m_sensorInd(1,s_measurementInd));
                v_rxPos = m_sensorPos(:,m_sensorInd(2,s_measurementInd));

                phi1 = norm(v_txPos-v_rxPos);
                phi2 = sqrt( (x_axis-v_txPos(1)).^2 + (y_axis-v_txPos(2)).^2 ) + sqrt( (x_axis-v_rxPos(1)).^2 + (y_axis-v_rxPos(2)).^2 );

                m_W = obj.h_w(phi1,phi2);
                m_vecWCollection(:,s_measurementInd) = m_W(:);
            end

            % 2. Estimate m_F_est according to obj.ch_reg_f_type (regularizer type) and ch_calibrationType.
            switch obj.ch_calibrationType
				case 'none'
					[v_sumOfGains,v_pathLoss] = obj.parameterEstimator(m_sensorPos,m_sensorInd,[]);
					s_check = -1 * (v_measurements - (v_sumOfGains - v_pathLoss));
					v_f_est = obj.chooseSolver(s_check,m_vecWCollection);
				case 'previous'
					[v_sumOfGains,v_pathLoss] = obj.parameterEstimator(m_sensorPos,m_sensorInd,v_measurementsNoShadowing);
					s_check = -1 * (v_measurements - (v_sumOfGains - v_pathLoss));
					v_f_est = obj.chooseSolver(s_check,m_vecWCollection);
				case 'simultaneous'
					% we need an estimator jointly estimating v_f_est,
					% pathloss exponent, and tx/rx gains only with shadow-faded channel
					% gain measurements
					s_sensorNum = size(m_sensorPos,2);
					v_sensorDistancesdB = zeros(s_measurementNum,1);
					
					m_Omega = obj.sensorMapOp(m_sensorInd,s_sensorNum);
					for s_measurementInd = 1: s_measurementNum
						v_txPos = m_sensorPos(:,m_sensorInd(1,s_measurementInd));
						v_rxPos = m_sensorPos(:,m_sensorInd(2,s_measurementInd));
						v_sensorDistancesdB(s_measurementInd) = 10 * log10(norm(v_txPos-v_rxPos));
					end
					
					m_barI = eye(s_measurementNum) - (1/(norm(v_sensorDistancesdB,2)^2))*(v_sensorDistancesdB)*(v_sensorDistancesdB');
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
		
		function [m_F_est,h_w_est] = blindEstimation(obj,m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing)
			% Estimate a spatial loss field in blind fashion with a
			% chosen regularizer for the SLF(tikhonov, l-1, and total variation).
			% In addition, a weight function is estimated with a kernel
			% method. In particular, the weight function lies in
			% reproducing kernel hilbert space. Complexity of the function
			% is controllable by a hilbert norm.
                        
            % 1. Compute features and correponding kernel matrix. For a
            % practical purpose, features are clustered and represented by
            % much less number of cluster centroids. Those
            % centroids are found by either 'random' sampling, or
            % well-known 'k-means'.
			
            [N_x, N_y] = size(obj.ini_F);
            s_gridNum = N_x * N_y;
            s_measurementNum = length(v_measurements);
			
			% check
			assert( obj.s_clusterNum <= s_gridNum*s_measurementNum );
			
			% Code here
			% obtain centroids with inputs of
			% m_sensorPos,m_sensorInd,obj.s_clusterNum, and
			% obj.ch_clustType.
			
			
			% Code here
			% find a kernel matrix with centroids obtained above and a
			% selected kernel function obj.h_kernel
			
			
			
            % 2. Estimate m_F_est according to obj.ch_reg_f_type (regularizer type) and ch_calibrationType
			%    and v_alpha.
				
			% alternating minimization setup
            prev_v_f = obj.ini_F(:);
            s_stoppingCriterion = 1e-5;
			s_iterationMax = 1e2;
            estimateDifference = inf;
			
            switch obj.ch_calibrationType
				case 'none'
					[v_sumOfGains,v_pathLoss] = obj.parameterEstimator(m_sensorPos,m_sensorInd,[]);
					s_check = -1 * (v_measurements - (v_sumOfGains - v_pathLoss));
					
					s_iterationNum = 1;
					while (estimateDifference > s_stoppingCriterion) && (s_iterationNum < s_iterationMax)
						
						% [S1] estimate v_alpha
						
						% [S2] estimate v_f_est
						estimateDifference = norm(prev_v_f - v_f_est,2);
						prev_v_f = v_f_est;
						s_iterationNum = s_iterationNum + 1;
					end
					
					
					v_f_est = obj.chooseSolver(s_check,m_vecWCollection);
				case 'previous'
					[v_sumOfGains,v_pathLoss] = obj.parameterEstimator(m_sensorPos,m_sensorInd,v_measurementsNoShadowing);
					s_check = -1 * (v_measurements - (v_sumOfGains - v_pathLoss));
					v_f_est = obj.chooseSolver(s_check,m_vecWCollection);
				case 'simultaneous'
					% we need an estimator jointly estimating v_f_est,
					% pathloss exponent, and tx/rx gains only with shadow-faded channel
					% gain measurements
					s_sensorNum = size(m_sensorPos,2);
					v_sensorDistancesdB = zeros(s_measurementNum,1);
					
					m_Omega = obj.sensorMapOp(m_sensorInd,s_sensorNum);
					for s_measurementInd = 1: s_measurementNum
						v_txPos = m_sensorPos(:,m_sensorInd(1,s_measurementInd));
						v_rxPos = m_sensorPos(:,m_sensorInd(2,s_measurementInd));
						v_sensorDistancesdB(s_measurementInd) = 10 * log10(norm(v_txPos-v_rxPos));
					end
					
					m_barI = eye(s_measurementNum) - (1/(norm(v_sensorDistancesdB,2)^2))*(v_sensorDistancesdB)*(v_sensorDistancesdB');
					m_tempOmega = m_barI * m_Omega;
					% add "1e-4 * eye(s_sensorNum)" to "m_tempOmega' * m_tempOmega" for matrix inversion
					m_barOmega = m_Omega * ((m_tempOmega' * m_tempOmega + 1e-4 * eye(s_sensorNum)) \ (m_tempOmega' * m_barI));
					m_tildeI = m_barI * (eye(s_measurementNum) - m_barOmega);
					v_barmeasurements = m_tildeI * v_measurements;
					m_barvecWCollection = -1 * m_vecWCollection * m_tildeI';
					v_f_est = obj.chooseSolver(v_barmeasurements,m_barvecWCollection);
			end
          
			% blind estimator outputs
			m_F_est = reshape(v_f_est,N_x,N_y);			
			h_w_est = @(input)  myRepThm(alpha,evl_pnt,input,myKfunc);
		end
		
		function [v_sumOfGains,v_pathLoss] = parameterEstimator(obj,m_sensorPos,m_sensorInd,v_measurementsNoShadowing)
			% Estimate(or assign) the sumOfGain and pathloss for every pair
			% of sensors when obj.ch_calibrationType is either "none", or
			% "previous".
			
			s_measurementNum = size(m_sensorInd,2);
			s_sensorNum = size(m_sensorPos,2);
			v_sensorDistancesdB = zeros(s_measurementNum,1);
			
			m_Omega = obj.sensorMapOp(m_sensorInd,s_sensorNum);
			for s_measurementInd = 1: s_measurementNum
				v_txPos = m_sensorPos(:,m_sensorInd(1,s_measurementInd));
				v_rxPos = m_sensorPos(:,m_sensorInd(2,s_measurementInd));
				v_sensorDistancesdB(s_measurementInd) = 10 * log10(norm(v_txPos-v_rxPos));
			end
			
			if  isempty(v_measurementsNoShadowing) % when ch_calibrationType = 'none'
				if isempty(obj.v_gains)
					v_gains_est = zeros(s_sensorNum,1); % by default
				else
					v_gains_est = obj.v_gains; 
				end
				s_pathLossExponent_est = obj.s_pathLossExponent;
			else % when ch_calibrationType = 'previous'
				% v_estGains / s_estPathLossExponent should be estimated
				m_regMat = [m_Omega,-1 * v_sensorDistancesdB];
				v_parameters = (m_regMat'*m_regMat + 1e-12* eye(s_sensorNum + 1))\(m_regMat'*v_measurementsNoShadowing);
				v_gains_est = v_parameters(1:s_sensorNum,1);
				s_pathLossExponent_est = v_parameters(s_sensorNum+1,1);
			end
						
			% v_estGains and s_estPathLossExponent should be known or
			% estimated before this step.			
			v_sumOfGains = m_Omega * v_gains_est;
			v_pathLoss = s_pathLossExponent_est * v_sensorDistancesdB;
		end
		
		function v_f_est = chooseSolver(obj,s_check,m_vecWCollection)
			% choose a solver among 'tikhonov'/'l1_ISTA'/'l1_PCO'/'totalvariation' 
			% according to a regularization type (ch_reg_f_type).
			[N_x, N_y] = size(obj.ini_F);
            switch obj.ch_reg_f_type

                case 'tikhonov'
                    v_f_est = ChannelGainMapEstimator.ridgeRegression(m_vecWCollection,s_check,obj.mu_f,obj.m_tikhonov);
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
            while (estimateDifference > s_stoppingCriterion) && (s_iterationNum < 1e4)

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
            v_seed2 = [1,zeros(1,N_x-1),-1,zeros(1,N_g - (N_x +1))];
            m_D2 = zeros(N_g-N_x,N_g);
            for i = 1 : N_g-N_x
                m_D2(i,:) = circshift(v_seed2,[1,i-1]);
            end

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
		
		function m_K = kernelMatrix()
		end

	end
	
	methods(Static) 
	
		function m_imageOut = postprocess(m_imageIn)
			
			m_imageOut = m_imageIn;
			m_imageOut( m_imageOut < 0 ) = 0;
			m_imageOut( m_imageOut > 1 ) = 1;
			
		end
		
	end
	
end

