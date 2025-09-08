classdef WavefrontAnalyzer < handle
    % Simple wavefront analysis class
    
    properties
        wavelength = 550e-9;  % meters
        grid_size = 256;
        wavefront_data
    end
    
    methods
        function obj = WavefrontAnalyzer(wavelength)
            % Constructor
            if nargin > 0
                obj.wavelength = wavelength;
            end
            fprintf('WavefrontAnalyzer initialized\n');
        end
        
        function generateTestData(obj)
            % Generate simple test wavefront
            [x, y] = meshgrid(linspace(-1,1,obj.grid_size));
            r = sqrt(x.^2 + y.^2);
            mask = r <= 1.0;
            
            % Simple defocus pattern
            obj.wavefront_data = 0.1 * (r.^2 - 0.5);
            obj.wavefront_data(~mask) = NaN;
        end
        
        function rms_error = calculateRMS(obj)
            % Calculate RMS wavefront error
            valid_data = obj.wavefront_data(~isnan(obj.wavefront_data));
            rms_error = sqrt(mean(valid_data.^2));
        end
        
        function plotWavefront(obj)
            % Simple wavefront plot
            imagesc(obj.wavefront_data);
            colorbar;
            title('Wavefront Map');
            axis equal tight;
        end
    end
end