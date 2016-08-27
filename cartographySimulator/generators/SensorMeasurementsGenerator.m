classdef SensorMeasurementsGenerator < Parameter
	% Subclasses 
	
	
	properties(Constant)
	end
	
	
		
	methods
		
		function obj = SensorMeasurementsGenerator(varargin) 
			obj@Parameter(varargin{:});
		end
		
	end
	

	
	methods(Abstract)
		
		[m_sensorPos,m_sensorInd,v_measurements] = realization(obj);
		% m_sensorPos  is a 2 x #sensors matrix where the i-th column
		%              indicates the (x,y) coordinates of the i-th sensor
		% m_sensorInd  is a 2 x #measurements matrix where m_sensorInd(1,i)
		%              and m_sensorInd(2,i) are the indices of the sensors
		%              for measurement i.
		% v_measurements is a 1 x #measurements vector where
		%              v_measurements(i) is the i-th measurement, containing
		%              path loss, the gain of sensors with indices m_sensorInd(1,i)
		%              and m_sensorInd(2,i), and the shadowing due to
		%              absorption of the medium. 
		%
		
		
	end
	
	
end

