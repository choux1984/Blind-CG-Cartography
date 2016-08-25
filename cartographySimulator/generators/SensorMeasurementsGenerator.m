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
		
		[s_check,Tx_pos,Rx_pos] = realization(obj);
		%   s_check       (obj.s_measurementNum)-by-1 vector of corrupted noise signals
		%   Tx_pos        2-by-( obj.s_measurementNum) collections of Tx coordinates
		%   Rx_pos        2-by-( obj.s_measurementNum) collections of Rx coordinates
		
	end
	
	
end

