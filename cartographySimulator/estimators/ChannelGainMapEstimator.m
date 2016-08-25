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
		             %     'l1' = l_1;
		             %     'totalvariation' = total variation				
		mu_f         % Regularization parameter for spatial loss field
				
		ch_estimationType = 'non-blind';
		             % 'non-blind'
					 % 'blind'
					 %
		
		% non-blind estimation
		h_w = [];    %  
		
		% blind estimation
		mu_w         % Regularization parameter for weight function		
		myKfunc = [];% Kernel function: must take a scalar ---> change to take two vectors
		Nc           % Number of clusters; Set Nc = [] for no clustering		
		lambda_W     % h_w is estimated with domain equal to an ellipse whose semi-minor axis is lambda_W
		rho          % step size for ADMM algorithm
	end
		
	methods
		
		function obj = ChannelGainMapEstimator(varargin)
			obj@Parameter(varargin{:});
		end
		
	
	end
	
	methods
				
		function [m_F_est,h_w_est,m_centroids] = estimate(obj,s_check,m_txPos,m_rxPos)	
			%
			% OUTPUT:
			%   m_est_F         Nx x Ny matrix with the estimate of F
			%   h_w_est    - non-blind estimation
			%                 h_w_est = []
			%              - blind estimation
			%                 h_w_est is the estimate of h_w_est
			%
			%  m_centroids     2-by-(n_measurements*Ng): collection of all feature vectors
			
			
			% dependent variables
			[N_x, N_y] = size(obj.ini_F);
			Ng = N_x * N_y; % # of grid points			
			t = length(s_check);
						
			% Estimation
			switch obj.ch_estimationType
				case 'non-blind'					
					m_F_est  = obj.nonblindEstimation(s_check,m_txPos,m_rxPos);										
					h_w_est = [];
					
				case 'blind'

					% (finish this)
					assert( obj.Nc <= Ng*t );
					[evl_pnt,idx_phi,phi_col]  = generate_feature_vecs(N_x,N_y,Tx_pos,Rx_pos,Nc,clustering_type,lambda_W);
					K = myKmatrix(evl_pnt,myKfunc);
					
					[alpha,est_f] = blind_est(K,idx_phi,s_check,mu_w,ini_F(:),eps_error,max_itr,mu_f,reg_f_type,rho);
					w = @(input)  myRepThm(alpha,evl_pnt,input,myKfunc);
					est_F = reshape(est_f,N_x,N_y);
				
				
					
				otherwise
					error('unrecognized option');
			end
				 
		end
		
		function m_F_est  = nonblindEstimation(obj,s_check,m_txPos,m_rxPos)
			% write this using obj.ch_reg_f_type, obj.h_w and obj.mu_f
			
			% check
			if isempty(obj.ini_F)&&(strcmp(obj.ch_reg_f_type,'tikhonov')==0)
				assert(~isempty(obj.ini_F));				
			end
						
			switch obj.ch_reg_f_type
				
				case 'tikhonov'
					m_F_est = [];
					
				case 'l1'
					m_F_est = [];
					
				case 'totalvariation'
					m_F_est = [];
			end
			
		end
		
	end
	
end

