classdef WavefrontAnalyzer < handle
    % WAVEFRONTANALYZER A MATLAB class for wavefront analysis and reconstruction
    %
    % This class provides methods for analyzing wavefront data from various
    % sources including Shack-Hartmann sensors, interferometers, and 
    % synthetic data generation for optical system modeling.
    %
    % Example usage:
    %   wfa = WavefrontAnalyzer();
    %   wfa.generateSyntheticWavefront('zernike', [0.1, 0.2, 0.05]);
    %   wfa.analyzeWavefront();
    %   wfa.plotResults();
    
    properties (Access = public)
        % Wavefront data and grid
        wavefront           % 2D wavefront map (waves)
        x_grid             % X coordinate grid (normalized)
        y_grid             % Y coordinate grid (normalized)
        aperture_mask      % Circular aperture mask
        
        % System parameters
        wavelength         % Wavelength in meters
        pupil_diameter     % Pupil diameter in meters
        grid_size          % Grid size (N x N)
        
        % Analysis results
        zernike_coeffs     % Zernike polynomial coefficients
        rms_error          % RMS wavefront error
        pv_error           % Peak-to-valley error
        strehl_ratio       % Strehl ratio
        
        % Plotting properties
        figure_handles     % Handles to created figures
    end
    
    properties (Access = private)
        % Internal calculation matrices
        zernike_matrix     % Matrix of Zernike polynomials
        radial_coords      % Radial coordinate array
        angular_coords     % Angular coordinate array
    end
    
    methods (Access = public)
        function obj = WavefrontAnalyzer(varargin)
            % Constructor for WavefrontAnalyzer
            % 
            % Optional name-value pairs:
            %   'wavelength' - wavelength in meters (default: 550e-9)
            %   'pupil_diameter' - pupil diameter in meters (default: 0.01)
            %   'grid_size' - grid size N for N x N array (default: 256)
            
            % Parse input arguments
            p = inputParser;
            addParameter(p, 'wavelength', 550e-9, @(x) x > 0);
            addParameter(p, 'pupil_diameter', 0.01, @(x) x > 0);
            addParameter(p, 'grid_size', 256, @(x) x > 0 && mod(x,1) == 0);
            parse(p, varargin{:});
            
            % Set properties from parsed inputs
            obj.wavelength = p.Results.wavelength;
            obj.pupil_diameter = p.Results.pupil_diameter;
            obj.grid_size = p.Results.grid_size;
            
            % Initialize coordinate grids
            obj.initializeGrids();
            
            % Initialize other properties
            obj.figure_handles = [];
            obj.zernike_coeffs = [];
            
            fprintf('WavefrontAnalyzer initialized:\n');
            fprintf('  Wavelength: %.0f nm\n', obj.wavelength * 1e9);
            fprintf('  Pupil diameter: %.1f mm\n', obj.pupil_diameter * 1000);
            fprintf('  Grid size: %d x %d\n', obj.grid_size, obj.grid_size);
        end
        
        function generateSyntheticWavefront(obj, method, parameters)
            % Generate synthetic wavefront data for testing
            %
            % Inputs:
            %   method - 'zernike', 'sine', or 'random'
            %   parameters - method-specific parameters
            
            switch lower(method)
                case 'zernike'
                    obj.generateZernikeWavefront(parameters);
                    
                case 'sine'
                    obj.generateSineWavefront(parameters);
                    
                case 'random'
                    obj.generateRandomWavefront(parameters);
                    
                otherwise
                    error('Unknown wavefront generation method: %s', method);
            end
            
            % Apply aperture mask
            obj.wavefront(~obj.aperture_mask) = NaN;
            
            fprintf('Synthetic wavefront generated using %s method\n', method);
        end
        
        function loadWavefrontData(obj, filename)
            % Load wavefront data from file
            %
            % Supported formats: .mat, .fits, .csv
            
            [~, ~, ext] = fileparts(filename);
            
            switch lower(ext)
                case '.mat'
                    data = load(filename);
                    if isfield(data, 'wavefront')
                        obj.wavefront = data.wavefront;
                    else
                        fields = fieldnames(data);
                        obj.wavefront = data.(fields{1});
                    end
                    
                case '.fits'
                    try
                        obj.wavefront = fitsread(filename);
                    catch
                        error('FITS file reading requires Image Processing Toolbox');
                    end
                    
                case '.csv'
                    obj.wavefront = readmatrix(filename);
                    
                otherwise
                    error('Unsupported file format: %s', ext);
            end
            
            % Resize if necessary
            if size(obj.wavefront, 1) ~= obj.grid_size
                obj.wavefront = imresize(obj.wavefront, [obj.grid_size, obj.grid_size]);
            end
            
            % Apply aperture mask
            obj.wavefront(~obj.aperture_mask) = NaN;
            
            fprintf('Wavefront data loaded from %s\n', filename);
        end
        
        function fitZernikePolynomials(obj, max_order)
            % Fit Zernike polynomials to the wavefront data
            %
            % Input:
            %   max_order - maximum radial order to fit (default: 6)
            
            if nargin < 2
                max_order = 6;
            end
            
            if isempty(obj.wavefront)
                error('No wavefront data available. Generate or load data first.');
            end
            
            % Generate Zernike polynomial matrix
            obj.generateZernikeMatrix(max_order);
            
            % Get valid data points (inside aperture)
            valid_mask = obj.aperture_mask & ~isnan(obj.wavefront);
            valid_wavefront = obj.wavefront(valid_mask);
            valid_zernikes = obj.zernike_matrix(valid_mask, :);
            
            % Least squares fit
            obj.zernike_coeffs = valid_zernikes \ valid_wavefront;
            
            fprintf('Zernike fit completed with %d coefficients\n', length(obj.zernike_coeffs));
            
            % Calculate fit quality
            fitted_wavefront = obj.zernike_matrix * obj.zernike_coeffs;
            fitted_wavefront = reshape(fitted_wavefront, size(obj.wavefront));
            residual = obj.wavefront - fitted_wavefront;
            
            residual_rms = sqrt(nanmean(residual(valid_mask).^2));
            original_rms = sqrt(nanmean(valid_wavefront.^2));
            
            fprintf('Fit quality: %.1f%% RMS reduction\n', ...
                    (1 - residual_rms/original_rms) * 100);
        end
        
        function results = analyzeWavefront(obj)
            % Perform comprehensive wavefront analysis
            %
            % Returns:
            %   results - structure containing analysis results
            
            if isempty(obj.wavefront)
                error('No wavefront data available. Generate or load data first.');
            end
            
            % Get valid data points
            valid_mask = obj.aperture_mask & ~isnan(obj.wavefront);
            valid_data = obj.wavefront(valid_mask);
            
            % Basic statistics
            obj.rms_error = sqrt(mean(valid_data.^2));
            obj.pv_error = max(valid_data) - min(valid_data);
            
            % Strehl ratio calculation (approximate)
            phase_variance = var(valid_data) * (2*pi)^2;
            obj.strehl_ratio = exp(-phase_variance);
            
            % Power spectral density
            [psd, freq] = obj.calculatePSD();
            
            % Gradient analysis
            [grad_x, grad_y] = gradient(obj.wavefront);
            grad_magnitude = sqrt(grad_x.^2 + grad_y.^2);
            max_gradient = nanmax(grad_magnitude(:));
            
            % Compile results
            results = struct();
            results.rms_error = obj.rms_error;
            results.pv_error = obj.pv_error;
            results.strehl_ratio = obj.strehl_ratio;
            results.max_gradient = max_gradient;
            results.psd = psd;
            results.frequency = freq;
            results.valid_points = sum(valid_mask(:));
            results.fill_factor = sum(valid_mask(:)) / numel(valid_mask);
            
            if ~isempty(obj.zernike_coeffs)
                results.zernike_coeffs = obj.zernike_coeffs;
                results.zernike_rms = sqrt(sum(obj.zernike_coeffs(2:end).^2));
            end
            
            fprintf('\nWavefront Analysis Results:\n');
            fprintf('  RMS Error: %.3f waves\n', results.rms_error);
            fprintf('  P-V Error: %.3f waves\n', results.pv_error);
            fprintf('  Strehl Ratio: %.3f\n', results.strehl_ratio);
            fprintf('  Fill Factor: %.1f%%\n', results.fill_factor * 100);
            
            if ~isempty(obj.zernike_coeffs)
                fprintf('  Zernike RMS: %.3f waves\n', results.zernike_rms);
            end
        end
        
        function plotResults(obj)
            % Create comprehensive plots of wavefront analysis results
            
            if isempty(obj.wavefront)
                error('No wavefront data available.');
            end
            
            % Main wavefront plot
            fig1 = figure('Name', 'Wavefront Analysis', 'NumberTitle', 'off');
            obj.figure_handles = [obj.figure_handles, fig1];
            
            % Wavefront map
            subplot(2, 3, 1);
            imagesc(obj.x_grid(1,:), obj.y_grid(:,1), obj.wavefront);
            axis equal tight;
            colorbar;
            title('Wavefront Map (waves)');
            xlabel('Normalized X');
            ylabel('Normalized Y');
            colormap(gca, 'RdBu_r');
            
            % 3D surface plot
            subplot(2, 3, 2);
            surf(obj.x_grid, obj.y_grid, obj.wavefront, 'EdgeColor', 'none');
            title('3D Wavefront Surface');
            xlabel('Normalized X');
            ylabel('Normalized Y');
            zlabel('Wavefront (waves)');
            colormap(gca, 'jet');
            
            % Cross-sections
            subplot(2, 3, 3);
            center = ceil(obj.grid_size / 2);
            x_profile = obj.wavefront(center, :);
            y_profile = obj.wavefront(:, center);
            plot(obj.x_grid(1,:), x_profile, 'r-', 'LineWidth', 2);
            hold on;
            plot(obj.y_grid(:,1), y_profile, 'b--', 'LineWidth', 2);
            title('Central Cross-sections');
            xlabel('Normalized Coordinate');
            ylabel('Wavefront (waves)');
            legend('X-profile', 'Y-profile', 'Location', 'best');
            grid on;
            
            % Histogram
            subplot(2, 3, 4);
            valid_data = obj.wavefront(obj.aperture_mask & ~isnan(obj.wavefront));
            histogram(valid_data, 50, 'Normalization', 'probability');
            title('Wavefront Distribution');
            xlabel('Wavefront (waves)');
            ylabel('Probability');
            grid on;
            
            % Power spectral density
            subplot(2, 3, 5);
            try
                [psd, freq] = obj.calculatePSD();
                loglog(freq, psd, 'b-', 'LineWidth', 1.5);
                title('Power Spectral Density');
                xlabel('Spatial Frequency (cycles/aperture)');
                ylabel('PSD (waves^2·cycles^{-1})');
                grid on;
            catch
                text(0.5, 0.5, 'PSD calculation failed', ...
                     'HorizontalAlignment', 'center');
                title('Power Spectral Density');
            end
            
            % Zernike coefficients (if available)
            subplot(2, 3, 6);
            if ~isempty(obj.zernike_coeffs)
                stem(1:length(obj.zernike_coeffs), obj.zernike_coeffs, 'filled');
                title('Zernike Coefficients');
                xlabel('Zernike Index');
                ylabel('Coefficient (waves RMS)');
                grid on;
            else
                text(0.5, 0.5, 'Fit Zernikes first', ...
                     'HorizontalAlignment', 'center');
                title('Zernike Coefficients');
            end
            
            % Adjust figure layout
            sgtitle(sprintf('Wavefront Analysis (λ = %.0f nm)', ...
                           obj.wavelength * 1e9), 'FontSize', 14, 'FontWeight', 'bold');
        end
        
        function exportResults(obj, filename, format)
            % Export analysis results to file
            %
            % Inputs:
            %   filename - output filename
            %   format - 'mat', 'csv', or 'excel'
            
            if nargin < 3
                format = 'mat';
            end
            
            % Prepare data structure
            export_data = struct();
            export_data.wavefront = obj.wavefront;
            export_data.x_grid = obj.x_grid;
            export_data.y_grid = obj.y_grid;
            export_data.wavelength = obj.wavelength;
            export_data.pupil_diameter = obj.pupil_diameter;
            
            if ~isempty(obj.zernike_coeffs)
                export_data.zernike_coeffs = obj.zernike_coeffs;
            end
            
            % Perform analysis if not done already
            if isempty(obj.rms_error)
                obj.analyzeWavefront();
            end
            
            export_data.rms_error = obj.rms_error;
            export_data.pv_error = obj.pv_error;
            export_data.strehl_ratio = obj.strehl_ratio;
            
            switch lower(format)
                case 'mat'
                    save(filename, '-struct', 'export_data');
                    
                case 'csv'
                    writematrix(obj.wavefront, filename);
                    
                case 'excel'
                    try
                        writetable(struct2table(export_data), filename);
                    catch
                        warning('Excel export failed. Saving as MAT file instead.');
                        save([filename, '.mat'], '-struct', 'export_data');
                    end
                    
                otherwise
                    error('Unsupported export format: %s', format);
            end
            
            fprintf('Results exported to %s\n', filename);
        end
    end
    
    methods (Access = private)
        function initializeGrids(obj)
            % Initialize coordinate grids and aperture mask
            
            [obj.x_grid, obj.y_grid] = meshgrid(...
                linspace(-1, 1, obj.grid_size), ...
                linspace(-1, 1, obj.grid_size));
            
            % Create circular aperture
            obj.aperture_mask = (obj.x_grid.^2 + obj.y_grid.^2) <= 1.0;
            
            % Convert to polar coordinates for Zernike calculations
            obj.radial_coords = sqrt(obj.x_grid.^2 + obj.y_grid.^2);
            obj.angular_coords = atan2(obj.y_grid, obj.x_grid);
        end
        
        function generateZernikeWavefront(obj, coefficients)
            % Generate wavefront from Zernike coefficients
            
            if isempty(coefficients)
                coefficients = [0, 0.1, 0.05, 0.2, 0.1, -0.05];
            end
            
            obj.generateZernikeMatrix(length(coefficients));
            obj.wavefront = reshape(obj.zernike_matrix * coefficients(:), ...
                                  [obj.grid_size, obj.grid_size]);
        end
        
        function generateSineWavefront(obj, parameters)
            % Generate sinusoidal wavefront pattern
            
            if isempty(parameters)
                parameters = [2, 3, 0.2]; % [freq_x, freq_y, amplitude]
            end
            
            freq_x = parameters(1);
            freq_y = parameters(2);
            amplitude = parameters(3);
            
            obj.wavefront = amplitude * sin(2*pi*freq_x*obj.x_grid) .* ...
                           sin(2*pi*freq_y*obj.y_grid);
        end
        
        function generateRandomWavefront(obj, parameters)
            % Generate random wavefront with specified statistics
            
            if isempty(parameters)
                parameters = [0.1, 0.5]; % [rms, correlation_length]
            end
            
            rms_amplitude = parameters(1);
            if length(parameters) > 1
                correlation_length = parameters(2);
            else
                correlation_length = 0.1;
            end
            
            % Generate correlated random field
            obj.wavefront = rms_amplitude * randn(obj.grid_size);
            
            % Apply spatial filtering for correlation
            sigma = correlation_length * obj.grid_size / 2;
            h = fspecial('gaussian', [obj.grid_size, obj.grid_size], sigma);
            obj.wavefront = ifft2(fft2(obj.wavefront) .* fft2(h, obj.grid_size, obj.grid_size));
            obj.wavefront = real(obj.wavefront);
            obj.wavefront = obj.wavefront * rms_amplitude / std(obj.wavefront(:));
        end
        
        function generateZernikeMatrix(obj, num_coeffs)
            % Generate matrix of Zernike polynomials
            
            % Zernike ordering (n, m pairs)
            zernike_orders = obj.getZernikeOrdering(num_coeffs);
            
            obj.zernike_matrix = zeros(numel(obj.x_grid), length(zernike_orders));
            
            for i = 1:length(zernike_orders)
                n = zernike_orders(i, 1);
                m = zernike_orders(i, 2);
                
                Z = obj.zernikePoly(n, m, obj.radial_coords, obj.angular_coords);
                obj.zernike_matrix(:, i) = Z(:);
            end
        end
        
        function orders = getZernikeOrdering(obj, num_coeffs)
            % Get Zernike (n, m) ordering up to specified number of coefficients
            
            orders = [];
            n = 0;
            count = 0;
            
            while count < num_coeffs
                for m = -n:2:n
                    if count >= num_coeffs
                        break;
                    end
                    orders = [orders; n, m];
                    count = count + 1;
                end
                n = n + 1;
            end
        end
        
        function Z = zernikePoly(obj, n, m, rho, theta)
            % Calculate Zernike polynomial
            
            % Radial polynomial
            R = zeros(size(rho));
            for k = 0:(n-abs(m))/2
                coeff = (-1)^k * factorial(n-k) / ...
                        (factorial(k) * factorial((n+abs(m))/2-k) * factorial((n-abs(m))/2-k));
                R = R + coeff * rho.^(n-2*k);
            end
            
            % Azimuthal component
            if m >= 0
                Z = R .* cos(m * theta);
            else
                Z = R .* sin(abs(m) * theta);
            end
            
            % Normalization
            if n > 0
                Z = Z * sqrt(2*(n+1));
            end
        end
        
        function [psd, freq] = calculatePSD(obj)
            % Calculate power spectral density of the wavefront
            
            % Remove NaNs by interpolation
            valid_mask = ~isnan(obj.wavefront) & obj.aperture_mask;
            if sum(valid_mask(:)) < 0.5 * numel(obj.wavefront)
                error('Too many invalid points for PSD calculation');
            end
            
            % Simple PSD calculation
            wf_fft = fft2(obj.wavefront - nanmean(obj.wavefront(:)));
            psd_2d = abs(wf_fft).^2;
            
            % Radial average
            [psd, freq] = obj.radialAverage(psd_2d);
        end
        
        function [profile, r] = radialAverage(obj, data)
            % Calculate radial average of 2D data
            
            center = ceil(size(data) / 2);
            [X, Y] = meshgrid(1:size(data,2), 1:size(data,1));
            R = sqrt((X - center(2)).^2 + (Y - center(1)).^2);
            
            max_r = min(center) - 1;
            r = 0:max_r;
            profile = zeros(size(r));
            
            for i = 1:length(r)
                mask = (R >= r(i)-0.5) & (R < r(i)+0.5);
                if any(mask(:))
                    profile(i) = nanmean(data(mask));
                end
            end
            
            % Convert to normalized frequency
            r = r / max_r;
        end
    end
end