classdef SyntheticSensorMeasurementsGenerator < SensorMeasurementsGenerator
	% 
	
	properties(Constant)
		
	end

	properties % Required by superclass Parameter
		c_parsToPrint    = {};
		c_stringToPrint  = {};
		c_patternToPrint = {};
	end
		
	
	properties
		m_F            % spatial loss field matrix
		h_w            % weight function handle; it is a function of two 
		               % parameters phi_1 and phi_2 (see writeup)		
		s_measurementNum   % number of measurements. Each measurement corresponds
		               % to two sensors selected among the length(v_gains)
		               % sensors
		v_gains = [];  % vector with the gains of the sensors. If empty, then 
		               % a pair of sensor locations is generated per
		               % measurement and all have gain 0dB. 
        s_pathLossExponent = 2;
		
		s_noiseVar
		
		ch_sensorLocationMode = 'uniform'; % other options are 'boundary'
			
	end
	
	methods
		
		function obj = SyntheticSensorMeasurementsGenerator(varargin) 
			obj@SensorMeasurementsGenerator(varargin{:});
		end
		
	end
	
	methods
				
		function [m_sensorPos,m_sensorInd,v_measurements] = realization(obj)
			% OUTPUT: see parent class
			
			% 1. Generate sensor locations
			if isempty(obj.v_gains)  % one pair of sensors generated per measurement
				
				switch(obj.ch_sensorLocationMode)
					case 'uniform'
						[N_x,N_y] = size(obj.m_F);
						m_txPos = diag([N_x,N_y])*rand(2,obj.s_measurementNum);
						m_rxPos = diag([N_x,N_y])*rand(2,obj.s_measurementNum);
						
					case 'boundary'
						% to be coded
				end
				
				% v_sumOfGains = zeros ...
			
			else  % the locations of length(obj.v_gains) sensors are generated and 
				  % obj.s_measurementNum different pairs of sensors are
				  % selected and a measurement is generated per pair. 
				
				
				 % v_sumOfGains = ...
			end
			
			% 2. Obtain measurements from obj.m_F, obj.h_w, and sensor
			% locations
			v_shadowing = zeros(size(m_txPos,2),1);
			for s_measurementInd = 1:size(m_txPos,2)
				v_shadowing(s_measurementInd) = SyntheticSensorMeasurementsGenerator.getNoiselessMeasurements(obj.m_F,obj.h_w,m_txPos(:,s_measurementInd),m_rxPos(:,s_measurementInd)) + sqrt(obj.s_noiseVar)*randn;
				%v_pathloss =  some function of m_txPos(:,s_measurementInd) and m_rxPos(:,s_measurementInd)
			end
			
			
			% code this
			% add path loss and gains
			v_measurements = v_pathLoss + v_sumOfGains + v_shadowing;
			
		end
		
		
	end
	
	methods(Static)
		
		function s = getNoiselessMeasurements(m_F,h_w,v_txPos,v_rxPos)
			%
			%
			% v_txPos :   2 x 1 vector with tx position
			% v_rxPos :   2 x 1 vector with rx position
			%
			% s :  shadowing measurement
			%
			
			% 1. compute m_W matrix
			[sizeX,sizeY] = size(m_F);
            x_axis = repmat((1:sizeY), sizeX, 1);    % x_axis of a grid
			y_axis = repmat((sizeX:-1:1)', 1, sizeY); % y_axis of a grid
            
            
			phi1 = norm(v_txPos-v_rxPos);
			phi2 = sqrt( (x_axis-v_txPos(1)).^2 + (y_axis-v_txPos(2)).^2 ) + sqrt( (x_axis-v_rxPos(1)).^2 + (y_axis-v_rxPos(2)).^2 );
			
%			% (non-efficient way; please optimize if needed)
% 			m_W = zeros(sizeX,sizeY);
% 			for c1 = 1:sizeX
% 				for c2 = 1:sizeY
% 					m_W(c1,c2) = h_w( phi1 , phi2(c1,c2) );
% 				end
%             end
			
            % Efficient way to calculate m_W
            m_W = h_w( phi1 , phi2);
						
			% 2. multiply m_F and m_W
			s = sum(sum(m_W.*m_F)); % signal w/o noise
						
			
	end
		
		
		
	end
	
end

