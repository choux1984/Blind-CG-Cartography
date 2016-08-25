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
		s_measurementNum      
		s_noiseVar
		
		ch_sensorLocationMode = 'uniform'; % other options are 'boudary'
			
	end
	
	methods
		
		function obj = SyntheticSensorMeasurementsGenerator(varargin) 
			obj@SensorMeasurementsGenerator(varargin{:});
		end
		
	end
	
	methods
				
		function [v_shadowingMeasuremens,m_txPos,m_rxPos] = realization(obj)
			% OUTPUT: see parent class
			
			
			% 1. Generate sensor locations
			switch(obj.ch_sensorLocationMode)
				case 'uniform'
					[N_x,N_y] = size(obj.m_F);
					m_txPos = diag([N_x-1,N_y-1])*rand(2,obj.s_measurementNum)+1;
					m_rxPos = diag([N_x-1,N_y-1])*rand(2,obj.s_measurementNum)+1;
					
				case 'boundary'
				    % to be coded			
			end
			
			
			% 2. Obtain measurements from obj.m_F, obj.h_w, and sensor
			% locations
			v_shadowingMeasuremens = zeros(size(m_txPos,2),1);
			for s_measurementInd = 1:size(m_txPos,2)
				v_shadowingMeasuremens(s_measurementInd) = SyntheticSensorMeasurementsGenerator.getNoiselessMeasurements(obj.m_F,obj.h_w,m_txPos(:,s_measurementInd),m_rxPos(:,s_measurementInd)) + sqrt(obj.s_noiseVar)*randn;
			end
			
			
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
			rows = repmat((1:sizeX)', 1, sizeY);  % first coordinate of grid point
			cols = repmat(1:sizeY, sizeX, 1);     % second coordinate of grid point
			
			phi1 = norm(v_txPos-v_rxPos);
			phi2 = sqrt( (rows-v_txPos(1)).^2 + (cols-v_txPos(2)).^2 ) + sqrt( (rows-v_rxPos(1)).^2 + (cols-v_rxPos(2)).^2 );
			
			% (non-efficient way; please optimize if needed)
			m_W = zeros(sizeX,sizeY);
			for c1 = 1:sizeX
				for c2 = 1:sizeY
					m_W(c1,c2) = h_w( phi1 , phi2(c1,c2) );
				end
			end
			
						
			% 2. multiply m_F and m_W
			s = sum(sum(m_W.*m_F)); % signal w/o noise
						
			
	end
		
		
		
	end
	
end

