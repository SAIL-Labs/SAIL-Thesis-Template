class OpticalSystem:
    """Simple optical system analysis class"""

    def __init__(self, focal_length=100.0):
        self.focal_length = focal_length
        self.aperture = 10.0
        self.wavelength = 0.55e-6

    def calculate_f_number(self):
        """Calculate f-number of the system"""
        return self.focal_length / self.aperture

    def diffraction_limit(self):
        """Calculate diffraction limited spot size"""
        return 1.22 * self.wavelength * self.calculate_f_number()

    def analyze(self):
        """Perform basic system analysis"""
        f_num = self.calculate_f_number()
        spot_size = self.diffraction_limit()

        print(f"F-number: {f_num:.1f}")
        print(f"Diffraction limit: {spot_size*1e6:.1f} um")

        return {"f_number": f_num, "spot_size": spot_size}


# Example usage
if __name__ == "__main__":
    system = OpticalSystem(focal_length=200.0)
    results = system.analyze()
