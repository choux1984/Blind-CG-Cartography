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
		v_gains        % vector with the gains of the sensors. If empty, then 
		               % a pair of sensor locations is generated per
		               % measurement and all have gain 0dB. 
        s_pathLossExponent
		
		s_noiseVar
		
		ch_sensorLocationMode = 'uniform'; % other options are 'boundary'
			
	end
	
	methods
		
		function obj = SyntheticSensorMeasurementsGenerator(varargin) 
			obj@SensorMeasurementsGenerator(varargin{:});
		end
		
	end
	
	methods % realization
				
		function [m_sensorPos,m_sensorInd,v_measurements,v_measurementsNoShadowing] = realization(obj)
			% OUTPUT: 
			%  m_sensorPos      (2-by-s_sensorNum) a matrix containing every position of a sensor
            %  m_sensorInd      (2-by-s_measurementNum) a matrix contating 
            %                    Tx/Rx sensor indices. The 1st row contains
            %                    indices of Transmitters. The 2nd row contains
            %                    indices of receivers.
            %  v_measurements    (s_measurementNum-by-1) a set of channel gain measurements.
            %  v_measurementsNoShadowing (s_measurementNum-by-1) a set of channel gain measurements from free space
            
            
			% 1. Generate sensor locations and sum of sensor gain on each
			% pair of sensors			
            [N_x,N_y] = size(obj.m_F);
			if isempty(obj.v_gains)  % one pair of sensors generated per measurement
				
				switch(obj.ch_sensorLocationMode)
					case 'uniform'						
						m_txPos = diag([N_x,N_y])*rand(2,obj.s_measurementNum);
						m_rxPos = diag([N_x,N_y])*rand(2,obj.s_measurementNum);
						m_sensorPos = [m_txPos,m_rxPos];
                        m_sensorInd = [1:obj.s_measurementNum; obj.s_measurementNum+1:2 * obj.s_measurementNum];
					case 'boundary'
						% to be coded
				end
				
				v_sumOfGains = zeros(obj.s_measurementNum,1);
			
			else  % the locations of length(obj.v_gains) sensors are generated and
				% obj.s_measurementNum different pairs of sensors are
				% selected and a measurement is generated per pair.
				switch(obj.ch_sensorLocationMode)
					case 'uniform'
						s_sensorNum = length(obj.v_gains);
						m_sensorPos = diag([N_x,N_y])*rand(2,s_sensorNum);
						
					case 'boundary'
						% to be coded
				end
				
				m_sensorInd = zeros(2,obj.s_measurementNum);
				v_sumOfGains = zeros(obj.s_measurementNum,1);
				for s_measurementInd = 1 : obj.s_measurementNum
					m_sensorInd(:,s_measurementInd) = randperm(s_sensorNum,2)';
					v_sumOfGains(s_measurementInd,1) = obj.v_gains(m_sensorInd(1,s_measurementInd)) + obj.v_gains(m_sensorInd(2,s_measurementInd));
				end
				
			end
			
			% 2. Obtain measurements from obj.m_F, obj.h_w, and sensor
			% locations
			v_shadowing = zeros(obj.s_measurementNum,1);
            v_pathLoss = zeros(obj.s_measurementNum,1);

			for s_measurementInd = 1: obj.s_measurementNum
                v_txPos = m_sensorPos(:,m_sensorInd(1,s_measurementInd));
                v_rxPos = m_sensorPos(:,m_sensorInd(2,s_measurementInd));
				v_shadowing(s_measurementInd) = SyntheticSensorMeasurementsGenerator.getNoiselessMeasurements(obj.m_F,obj.h_w,v_txPos,v_rxPos);
                v_pathLoss(s_measurementInd) = 10 * obj.s_pathLossExponent * log10(norm(v_txPos-v_rxPos));
			end
			
			% add path loss and gains
			v_measurementsNoShadowing = v_sumOfGains - v_pathLoss + sqrt(obj.s_noiseVar)*randn;% free space measurements for calibration
			v_measurements =  v_sumOfGains - v_pathLoss  - v_shadowing + sqrt(obj.s_noiseVar)*randn; % put - sign on v_shadowing since the ground truth SLF only has nonnegative values.
			
		end
			
    end 
	
	methods(Static)
		
		function s_shadowMeasurement = getNoiselessMeasurements(m_F,h_w,v_txPos,v_rxPos)
			%
			% INPUT:
			%  v_txPos                  2 x 1 vector with tx position
			%  v_rxPos                  2 x 1 vector with rx position
			%
            % OUTPUT:
			%  s_shadowMeasurement      shadowing measurement
			%
			
			% 1. compute m_W matrix
			[N_x,N_y] = size(m_F);
            x_axis = repmat((0:N_y-1), N_x, 1);    % x_axis of a grid
			y_axis = repmat((N_x-1:-1:0)', 1, N_y); % y_axis of a grid
			           
			phi1 = norm(v_txPos-v_rxPos);
			phi2 = sqrt( (x_axis-v_txPos(1)).^2 + (y_axis-v_txPos(2)).^2 ) + sqrt( (x_axis-v_rxPos(1)).^2 + (y_axis-v_rxPos(2)).^2 );
						
            % Efficient way to calculate m_W
            m_W = h_w( phi1 , phi2);
			% 2. multiply m_F and m_W
			s_shadowMeasurement = sum(sum(m_W.*m_F)); % signal w/o noise
									
		end
						
	end
	
end

