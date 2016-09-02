classdef RealSensorMeasurementsGenerator < SensorMeasurementsGenerator
	% 
	
	properties(Constant)
		
	end

	properties % Required by superclass Parameter
		c_parsToPrint    = {};
		c_stringToPrint  = {};
		c_patternToPrint = {};
	end
			
	methods
		
		function obj = RealSensorMeasurementsGenerator(varargin) 
			obj@SensorMeasurementsGenerator(varargin{:});
		end
		
	end
	
	methods % realization
				
		function [t_sensorPos,t_sensorInd,v_measurements,v_measurementsNoShadowing] = realization(obj)
			% OUTPUT: 
			%  t_sensorPos      (2-by-s_sensorNum)-by-2 a tensor containing
			%                    every position of a sensor. The 1st slab
			%                    is for the sensors used in non-structured
			%                    test. The 2nd slab is for those used in
			%                    structured test.
            %  t_sensorInd      (2-by-s_measurementNum)-by-2 a tensor contating 
            %                    Tx/Rx sensor indices. The 1st row contains
            %                    indices of Transmitters. The 2nd row contains
            %                    indices of receivers. The 1st slab
			%                    is for the sensor indices used in non-structured
			%                    test. The 2nd slab is for those used in
			%                    structured test.
            %  v_measurements    (s_measurementNum-by-1) a set of channel gain measurements.
            %  v_measurementsNoShadowing (s_measurementNum-by-1) a set of channel gain measurements from free space
            
            
			% 1. Load sensor locations and indices for every pair of
			% sensors
			
			s_timeslot = 24;
			s_sensorNumPerSlot = 10;
			s_measurementNumPerSlot = 100;
			
			m_txPosNonshadowing = csvread('Tx_position_free.csv');
			m_rxPosNonshadowing = csvread('Rx_position_free.csv');
			m_txPos = csvread('Tx_position.csv');
			m_rxPos = csvread('Rx_position.csv');
			
			t_sensorPos(:,:,1) = [m_txPosNonshadowing',m_rxPosNonshadowing'];
			t_sensorPos(:,:,2) = [m_txPos',m_rxPos'];
			
			for s_timeIdx = 1 : s_timeslot		
				v_tempTxIdx = 1 + (s_timeIdx-1)*s_sensorNumPerSlot: s_timeIdx*s_sensorNumPerSlot;
				v_tempRxIdx = (s_timeslot*s_sensorNumPerSlot + 1) + (s_timeIdx-1)*s_sensorNumPerSlot: s_timeslot*s_sensorNumPerSlot + s_timeIdx*s_sensorNumPerSlot;

				t_sensorInd(1,1 + (s_timeIdx-1)*s_measurementNumPerSlot:s_timeIdx*s_measurementNumPerSlot) = kron(ones(1,s_sensorNumPerSlot),v_tempTxIdx);
				t_sensorInd(2,1 + (s_timeIdx-1)*s_measurementNumPerSlot:s_timeIdx*s_measurementNumPerSlot) = kron(v_tempRxIdx,ones(1,s_sensorNumPerSlot));
			end
			t_sensorInd(:,:,2) = t_sensorInd;
			
			% 2. Obtain measurements from obj.m_F, obj.h_w, and sensor
			% locations
			
			v_measurementsNoShadowing = 10 * log10(csvread('total_data_free.csv'));
			v_measurements = 10 * log10(csvread('total_data.csv'));
			
		end
		
		
    end 
	
end