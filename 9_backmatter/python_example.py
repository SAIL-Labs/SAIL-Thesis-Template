import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

class OpticalSystem:
    """
    A class for analyzing optical systems and wavefront aberrations.
    
    This class provides methods for calculating Zernike polynomials,
    analyzing point spread functions, and optimizing optical parameters.
    """
    
    def __init__(self, aperture_diameter=10.0, focal_length=100.0):
        """
        Initialize the optical system.
        
        Parameters:
        aperture_diameter (float): Diameter of the entrance pupil in mm
        focal_length (float): Focal length of the system in mm
        """
        self.aperture_diameter = aperture_diameter
        self.focal_length = focal_length
        self.f_number = focal_length / aperture_diameter
        self.wavelength = 0.55e-6  # Default wavelength in meters
        self.zernike_coefficients = np.zeros(15)  # Up to 4th order
        
    def calculate_zernike_polynomial(self, n, m, rho, theta):
        """
        Calculate Zernike polynomial Z_n^m at given polar coordinates.
        
        Parameters:
        n (int): Radial order
        m (int): Azimuthal order
        rho (array): Normalized radial coordinates (0 to 1)
        theta (array): Azimuthal angles in radians
        
        Returns:
        array: Zernike polynomial values
        """
        # Radial polynomial calculation
        R = np.zeros_like(rho)
        for k in range((n - abs(m)) // 2 + 1):
            coeff = (-1)**k * np.math.factorial(n - k)
            coeff /= (np.math.factorial(k) * 
                     np.math.factorial((n + abs(m)) // 2 - k) * 
                     np.math.factorial((n - abs(m)) // 2 - k))
            R += coeff * rho**(n - 2*k)
        
        # Apply azimuthal component
        if m >= 0:
            Z = R * np.cos(m * theta)
        else:
            Z = R * np.sin(abs(m) * theta)
            
        return Z
    
    def compute_wavefront(self, grid_size=256):
        """
        Compute the wavefront map from Zernike coefficients.
        
        Parameters:
        grid_size (int): Size of the square grid
        
        Returns:
        tuple: (x, y, wavefront) coordinate arrays and wavefront map
        """
        # Create coordinate grids
        x = np.linspace(-1, 1, grid_size)
        y = np.linspace(-1, 1, grid_size)
        X, Y = np.meshgrid(x, y)
        
        # Convert to polar coordinates
        rho = np.sqrt(X**2 + Y**2)
        theta = np.arctan2(Y, X)
        
        # Create circular aperture mask
        mask = rho <= 1.0
        
        # Calculate wavefront from Zernike polynomials
        wavefront = np.zeros_like(rho)
        zernike_orders = [(0, 0), (1, 1), (1, -1), (2, 0), (2, -2), (2, 2),
                         (3, -1), (3, 1), (3, -3), (3, 3), (4, 0), (4, 2),
                         (4, -2), (4, 4), (4, -4)]
        
        for i, (n, m) in enumerate(zernike_orders):
            if i < len(self.zernike_coefficients):
                Z = self.calculate_zernike_polynomial(n, m, rho, theta)
                wavefront += self.zernike_coefficients[i] * Z
        
        # Apply aperture mask
        wavefront[~mask] = np.nan
        
        return X, Y, wavefront
    
    def calculate_psf(self, grid_size=256, oversampling=4):
        """
        Calculate the point spread function from the wavefront.
        
        Parameters:
        grid_size (int): Size of the pupil grid
        oversampling (int): Oversampling factor for PSF calculation
        
        Returns:
        tuple: (u, v, psf) coordinate arrays and PSF
        """
        # Get wavefront
        X, Y, wavefront = self.compute_wavefront(grid_size)
        
        # Create complex amplitude in pupil plane
        amplitude = np.ones_like(wavefront)
        amplitude[np.isnan(wavefront)] = 0
        wavefront[np.isnan(wavefront)] = 0
        
        pupil = amplitude * np.exp(1j * 2 * np.pi * wavefront / self.wavelength)
        
        # Zero-pad for oversampling
        padded_size = grid_size * oversampling
        padded_pupil = np.zeros((padded_size, padded_size), dtype=complex)
        start = (padded_size - grid_size) // 2
        end = start + grid_size
        padded_pupil[start:end, start:end] = pupil
        
        # Calculate PSF via FFT
        psf_complex = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(padded_pupil)))
        psf = np.abs(psf_complex)**2
        
        # Create coordinate arrays (angular coordinates)
        angular_scale = self.wavelength / (2 * self.aperture_diameter / 2e3)  # radians
        u = np.linspace(-padded_size//2, padded_size//2-1, padded_size) * angular_scale
        v = u.copy()
        
        return u, v, psf
    
    def optimize_focus(self, target_strehl=0.8):
        """
        Optimize focus to maximize Strehl ratio.
        
        Parameters:
        target_strehl (float): Target Strehl ratio
        
        Returns:
        float: Optimal focus position (defocus Zernike coefficient)
        """
        def objective(defocus):
            self.zernike_coefficients[3] = defocus[0]  # Z_2^0 (defocus)
            _, _, psf = self.calculate_psf()
            strehl = np.max(psf) / np.sum(psf)  # Approximation
            return -strehl  # Minimize negative Strehl
        
        result = minimize(objective, [0.0], method='Nelder-Mead')
        optimal_defocus = result.x[0]
        self.zernike_coefficients[3] = optimal_defocus
        
        return optimal_defocus
    
    def set_aberrations(self, defocus=0, astigmatism_0=0, astigmatism_45=0,
                       coma_x=0, coma_y=0, spherical=0):
        """
        Set common aberration coefficients.
        
        Parameters:
        defocus (float): Defocus aberration (waves RMS)
        astigmatism_0 (float): 0-degree astigmatism (waves RMS)
        astigmatism_45 (float): 45-degree astigmatism (waves RMS)
        coma_x (float): Coma in x-direction (waves RMS)
        coma_y (float): Coma in y-direction (waves RMS)
        spherical (float): Spherical aberration (waves RMS)
        """
        self.zernike_coefficients[3] = defocus      # Z_2^0
        self.zernike_coefficients[4] = astigmatism_45  # Z_2^{-2}
        self.zernike_coefficients[5] = astigmatism_0   # Z_2^2
        self.zernike_coefficients[6] = coma_y       # Z_3^{-1}
        self.zernike_coefficients[7] = coma_x       # Z_3^1
        self.zernike_coefficients[10] = spherical   # Z_4^0
    
    def analyze_performance(self):
        """
        Analyze the optical system performance.
        
        Returns:
        dict: Dictionary containing performance metrics
        """
        # Calculate wavefront statistics
        X, Y, wavefront = self.compute_wavefront()
        valid_wavefront = wavefront[~np.isnan(wavefront)]
        
        # Calculate PSF
        u, v, psf = self.calculate_psf()
        
        # Performance metrics
        results = {
            'rms_wavefront_error': np.std(valid_wavefront),
            'pv_wavefront_error': np.max(valid_wavefront) - np.min(valid_wavefront),
            'strehl_ratio': np.max(psf) / np.sum(psf),
            'diffraction_limit': 1.22 * self.wavelength / self.aperture_diameter * 1e6,  # arcsec
            'f_number': self.f_number,
            'focal_length': self.focal_length,
            'aperture_diameter': self.aperture_diameter
        }
        
        return results

# Example usage and demonstration
if __name__ == "__main__":
    # Create an optical system
    system = OpticalSystem(aperture_diameter=200.0, focal_length=2000.0)
    
    # Add some aberrations
    system.set_aberrations(defocus=0.1, astigmatism_0=0.05, spherical=0.02)
    
    # Analyze performance
    performance = system.analyze_performance()
    print("Optical System Performance:")
    for key, value in performance.items():
        print(f"{key}: {value:.4f}")
    
    # Optimize focus
    optimal_focus = system.optimize_focus()
    print(f"\nOptimal focus adjustment: {optimal_focus:.4f} waves")
    
    # Generate plots (if running interactively)
    try:
        # Plot wavefront
        X, Y, wavefront = system.compute_wavefront()
        plt.figure(figsize=(12, 4))
        
        plt.subplot(131)
        plt.contourf(X, Y, wavefront, levels=20, cmap='RdBu_r')
        plt.colorbar(label='Wavefront Error (waves)')
        plt.title('Wavefront Map')
        plt.axis('equal')
        
        # Plot PSF
        u, v, psf = system.calculate_psf()
        plt.subplot(132)
        plt.imshow(psf, extent=[u[0], u[-1], v[0], v[-1]], cmap='hot')
        plt.colorbar(label='Intensity')
        plt.title('Point Spread Function')
        plt.xlabel('Angular coordinate (rad)')
        plt.ylabel('Angular coordinate (rad)')
        
        # Plot PSF cross-section
        plt.subplot(133)
        center = len(u) // 2
        plt.plot(u, psf[center, :] / np.max(psf))
        plt.xlabel('Angular coordinate (rad)')
        plt.ylabel('Normalized Intensity')
        plt.title('PSF Cross-section')
        plt.grid(True)
        
        plt.tight_layout()
        plt.show()
        
    except ImportError:
        print("Matplotlib not available for plotting")