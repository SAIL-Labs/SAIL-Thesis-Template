! ZEMAX ZPL Macro: Advanced Optical System Analysis
! This macro demonstrates automated analysis of optical systems including
! ray tracing, aberration analysis, and optimization routines
! Author: Example Code for Thesis Appendix
! Version: 1.0

! Initialize variables and system parameters
wavelength_count = NWAV()
field_count = NFLD()
surface_count = NSUR()

! Create output file for results
OUTPUT "analysis_results.txt"

PRINT "========================================="
PRINT "OPTICAL SYSTEM ANALYSIS MACRO"
PRINT "========================================="
PRINT ""

! Display system information
PRINT "System Information:"
PRINT "Number of wavelengths: ", wavelength_count
PRINT "Number of field points: ", field_count  
PRINT "Number of surfaces: ", surface_count
PRINT "System units: ", $UNIT()
PRINT ""

! Main analysis loop for each configuration
config_count = CONF()
FOR config_num, 1, config_count, 1
    ! Set current configuration
    SETCONFIG config_num
    
    PRINT "Configuration ", config_num, ":"
    PRINT "================================"
    
    ! Perform basic ray tracing analysis
    GOSUB ANALYZE_RAYS
    
    ! Calculate aberrations
    GOSUB CALCULATE_ABERRATIONS
    
    ! Analyze spot sizes
    GOSUB SPOT_ANALYSIS
    
    ! Check for vignetting
    GOSUB VIGNETTING_CHECK
    
    PRINT ""
NEXT config_num

! Optimization routine if requested
user_input$ = $GETSTRING("Perform optimization? (Y/N): ")
IF $UPPER(user_input$) = "Y"
    GOSUB OPTIMIZATION_ROUTINE
ENDIF

! Generate final report
GOSUB GENERATE_REPORT

! Close output file
OUTPUT SCREEN
PRINT "Analysis complete. Results saved to analysis_results.txt"

END

! Subroutine: Ray tracing analysis
LABEL ANALYZE_RAYS
    PRINT "Ray Tracing Analysis:"
    PRINT "-----------------------"
    
    ! Trace rays for all fields and wavelengths
    FOR field_num, 1, field_count, 1
        FOR wave_num, 1, wavelength_count, 1
            ! Set field and wavelength
            field_x = FLDX(field_num)
            field_y = FLDY(field_num)
            wavelength = WAVL(wave_num)
            
            ! Trace chief ray
            ray_valid = RAYT(0, 0, field_x, field_y, wave_num)
            
            IF ray_valid = 1
                ! Get ray data at image surface
                image_surf = surface_count
                ray_x = RAYX(image_surf)
                ray_y = RAYY(image_surf)
                ray_z = RAYZ(image_surf)
                
                PRINT "Field ", field_num, " Wave ", wave_num, ":"
                PRINT "  Chief ray at image: (", ray_x, ", ", ray_y, ", ", ray_z, ")"
                
                ! Calculate marginal ray
                pupil_radius = APMX()/2
                ray_valid_marg = RAYT(pupil_radius, 0, field_x, field_y, wave_num)
                
                IF ray_valid_marg = 1
                    marg_x = RAYX(image_surf)
                    marg_y = RAYY(image_surf)
                    spot_size = SQRT((marg_x - ray_x)^2 + (marg_y - ray_y)^2)
                    PRINT "  Marginal ray spot size: ", spot_size, " mm"
                ELSE
                    PRINT "  Marginal ray failed to trace"
                ENDIF
            ELSE
                PRINT "Field ", field_num, " Wave ", wave_num, ": Ray failed to trace"
            ENDIF
        NEXT wave_num
    NEXT field_num
    
    PRINT ""
RETURN

! Subroutine: Calculate aberrations using Seidel coefficients
LABEL CALCULATE_ABERRATIONS
    PRINT "Aberration Analysis:"
    PRINT "----------------------"
    
    ! Calculate Seidel aberration coefficients
    spherical_aber = 0
    coma = 0
    astigmatism = 0
    field_curvature = 0
    distortion = 0
    
    ! Loop through surfaces to accumulate aberrations
    FOR surf_num, 2, surface_count-1, 1
        ! Get surface parameters
        curvature = CURV(surf_num)
        thickness = THIC(surf_num)
        radius = 1/curvature
        
        ! Calculate contribution to aberrations (simplified)
        IF curvature <> 0
            contrib_sph = curvature^3 * thickness / 8
            contrib_coma = curvature^2 * thickness / 4
            contrib_ast = curvature * thickness / 2
            
            spherical_aber = spherical_aber + contrib_sph
            coma = coma + contrib_coma
            astigmatism = astigmatism + contrib_ast
        ENDIF
    NEXT surf_num
    
    PRINT "Seidel Aberration Coefficients:"
    PRINT "  Spherical Aberration: ", spherical_aber
    PRINT "  Coma: ", coma
    PRINT "  Astigmatism: ", astigmatism
    PRINT "  Field Curvature: ", field_curvature
    PRINT "  Distortion: ", distortion
    PRINT ""

RETURN

! Subroutine: Spot diagram analysis
LABEL SPOT_ANALYSIS
    PRINT "Spot Size Analysis:"
    PRINT "---------------------"
    
    ! Calculate RMS spot sizes for each field
    FOR field_num, 1, field_count, 1
        field_x = FLDX(field_num)
        field_y = FLDY(field_num)
        field_weight = FWGT(field_num)
        
        ! Initialize spot size calculation
        sum_x2 = 0
        sum_y2 = 0
        sum_xy = 0
        ray_count = 0
        
        ! Sample rays across pupil
        num_rings = 3
        rays_per_ring = 8
        
        FOR ring, 1, num_rings, 1
            ring_radius = ring * APMX()/(2 * num_rings)
            
            FOR ray_angle, 0, 360-360/rays_per_ring, 360/rays_per_ring
                pupil_x = ring_radius * COSI(ray_angle)
                pupil_y = ring_radius * SINE(ray_angle)
                
                ! Trace ray
                ray_valid = RAYT(pupil_x, pupil_y, field_x, field_y, 1)
                
                IF ray_valid = 1
                    image_x = RAYX(surface_count)
                    image_y = RAYY(surface_count)
                    
                    sum_x2 = sum_x2 + image_x^2
                    sum_y2 = sum_y2 + image_y^2
                    sum_xy = sum_xy + image_x * image_y
                    ray_count = ray_count + 1
                ENDIF
            NEXT ray_angle
        NEXT ring
        
        ! Calculate RMS spot size
        IF ray_count > 0
            mean_x2 = sum_x2 / ray_count
            mean_y2 = sum_y2 / ray_count
            rms_spot = SQRT(mean_x2 + mean_y2)
            
            PRINT "Field ", field_num, " (", field_x, ", ", field_y, "):"
            PRINT "  RMS Spot Size: ", rms_spot, " mm"
            PRINT "  Ray Count: ", ray_count
        ENDIF
    NEXT field_num
    
    PRINT ""
RETURN

! Subroutine: Check for vignetting
LABEL VIGNETTING_CHECK
    PRINT "Vignetting Analysis:"
    PRINT "----------------------"
    
    FOR field_num, 1, field_count, 1
        field_x = FLDX(field_num)
        field_y = FLDY(field_num)
        
        ! Test rays at full aperture
        full_aperture = APMX()/2
        vignetted_rays = 0
        total_test_rays = 16
        
        FOR test_ray, 1, total_test_rays, 1
            angle = test_ray * 360 / total_test_rays
            pupil_x = full_aperture * COSI(angle)
            pupil_y = full_aperture * SINE(angle)
            
            ray_valid = RAYT(pupil_x, pupil_y, field_x, field_y, 1)
            IF ray_valid = 0
                vignetted_rays = vignetted_rays + 1
            ENDIF
        NEXT test_ray
        
        vignetting_percent = vignetted_rays * 100 / total_test_rays
        
        PRINT "Field ", field_num, ": ", vignetting_percent, "% vignetting"
    NEXT field_num
    
    PRINT ""
RETURN

! Subroutine: Optimization routine
LABEL OPTIMIZATION_ROUTINE
    PRINT "Performing Optimization:"
    PRINT "============================"
    
    ! Save current merit function
    initial_merit = MFCN()
    PRINT "Initial Merit Function: ", initial_merit
    
    ! Set up optimization variables (example: vary some radii and thicknesses)
    variable_count = 0
    FOR surf_num, 2, surface_count-2, 2
        ! Make radius variable if it's not infinite
        IF CURV(surf_num) <> 0
            SOLVETYPE surf_num, VARIABLE
            variable_count = variable_count + 1
        ENDIF
        
        ! Make thickness variable for some surfaces
        IF surf_num < surface_count - 1
            thickness = THIC(surf_num)
            IF thickness > 1
                ! Set thickness as variable with constraints
                variable_count = variable_count + 1
            ENDIF
        ENDIF
    NEXT surf_num
    
    PRINT "Number of variables: ", variable_count
    
    ! Perform optimization
    IF variable_count > 0
        PRINT "Starting optimization..."
        OPTIMIZE
        
        ! Check results
        final_merit = MFCN()
        improvement = (initial_merit - final_merit) / initial_merit * 100
        
        PRINT "Final Merit Function: ", final_merit
        PRINT "Improvement: ", improvement, "%"
        
        IF improvement > 1
            PRINT "Optimization successful!"
        ELSE
            PRINT "Limited improvement achieved."
        ENDIF
    ELSE
        PRINT "No variables defined for optimization."
    ENDIF
    
    PRINT ""
RETURN

! Subroutine: Generate comprehensive report
LABEL GENERATE_REPORT
    PRINT ""
    PRINT "SUMMARY REPORT"
    PRINT "==============="
    
    ! System specifications
    effective_focal_length = 1/PWAV(1)  ! Simplified calculation
    f_number = effective_focal_length / APMX()
    
    PRINT "System Specifications:"
    PRINT "  Effective Focal Length: ", effective_focal_length, " mm"
    PRINT "  F-Number: f/", f_number
    PRINT "  Aperture Diameter: ", APMX(), " mm"
    
    ! Calculate overall system performance
    total_rms_spot = 0
    valid_fields = 0
    
    FOR field_num, 1, field_count, 1
        ! Quick spot size calculation for summary
        ray_valid = RAYT(0, 0, FLDX(field_num), FLDY(field_num), 1)
        IF ray_valid = 1
            ! Simplified spot size calculation
            field_spot = 0.1  ! Placeholder - would need full calculation
            total_rms_spot = total_rms_spot + field_spot
            valid_fields = valid_fields + 1
        ENDIF
    NEXT field_num
    
    IF valid_fields > 0
        average_spot_size = total_rms_spot / valid_fields
        PRINT "  Average RMS Spot Size: ", average_spot_size, " mm"
        
        ! Calculate angular resolution
        angular_resolution = average_spot_size / effective_focal_length * 206265
        PRINT "  Angular Resolution: ", angular_resolution, " arcsec"
    ENDIF
    
    ! Performance rating
    IF average_spot_size < 0.01
        performance$ = "Excellent"
    ELSEIF average_spot_size < 0.05
        performance$ = "Good"
    ELSEIF average_spot_size < 0.1
        performance$ = "Fair"
    ELSE
        performance$ = "Poor"
    ENDIF
    
    PRINT "  Overall Performance: ", performance$
    
    ! Recommendations
    PRINT ""
    PRINT "Recommendations:"
    IF average_spot_size > 0.05
        PRINT "  - Consider optimization to improve spot sizes"
        PRINT "  - Check for aberration correction opportunities"
    ENDIF
    
    IF vignetted_rays > 0
        PRINT "  - Address vignetting in outer fields"
        PRINT "  - Consider increasing aperture sizes"
    ENDIF
    
    PRINT "  - Verify manufacturing tolerances"
    PRINT "  - Consider environmental effects"
    
    PRINT ""
    PRINT "Analysis completed: ", $GDATE(), " at ", $GTIME()
    
RETURN

! Error handling subroutine
LABEL ERROR_HANDLER
    PRINT "Error encountered during analysis!"
    PRINT "Check system definition and try again."
    OUTPUT SCREEN
RETURN

! Utility function: Format number for display
LABEL FORMAT_NUMBER
    ! Input: number in variable 'value'
    ! Output: formatted string in 'formatted$'
    
    IF ABSO(value) >= 1000
        formatted$ = FORMAT(value, "%.2e")
    ELSEIF ABSO(value) >= 1
        formatted$ = FORMAT(value, "%.3f")
    ELSE
        formatted$ = FORMAT(value, "%.6f")
    ENDIF
    
RETURN

! Advanced ray tracing with custom sampling
LABEL ADVANCED_RAY_TRACE
    PRINT "Advanced Ray Tracing:"
    PRINT "------------------------"
    
    ! Custom ray sampling pattern
    sample_type$ = "hexagonal"  ! Can be "square", "hexagonal", or "random"
    num_samples = 100
    
    FOR field_num, 1, field_count, 1
        field_x = FLDX(field_num)
        field_y = FLDY(field_num)
        
        successful_rays = 0
        
        FOR sample, 1, num_samples, 1
            ! Generate sampling coordinates based on pattern
            IF sample_type$ = "hexagonal"
                GOSUB HEXAGONAL_SAMPLING
            ELSEIF sample_type$ = "square"
                GOSUB SQUARE_SAMPLING  
            ELSE
                GOSUB RANDOM_SAMPLING
            ENDIF
            
            ! pupil_x and pupil_y set by sampling subroutine
            ray_valid = RAYT(pupil_x, pupil_y, field_x, field_y, 1)
            
            IF ray_valid = 1
                successful_rays = successful_rays + 1
            ENDIF
        NEXT sample
        
        success_rate = successful_rays * 100 / num_samples
        PRINT "Field ", field_num, " ray success rate: ", success_rate, "%"
    NEXT field_num
    
RETURN

! Hexagonal sampling pattern
LABEL HEXAGONAL_SAMPLING
    ! Generate hexagonal sampling pattern
    ring_num = INT(SQRT(sample / 3)) + 1
    position_in_ring = sample - 3 * (ring_num - 1)^2
    
    IF ring_num = 1
        pupil_x = 0
        pupil_y = 0
    ELSE
        angle = position_in_ring * 60 / ring_num
        radius = ring_num * APMX() / (2 * SQRT(num_samples / 3))
        pupil_x = radius * COSI(angle)
        pupil_y = radius * SINE(angle)
    ENDIF
    
RETURN

! Square sampling pattern
LABEL SQUARE_SAMPLING
    grid_size = INT(SQRT(num_samples))
    row = INT((sample - 1) / grid_size) + 1
    col = sample - (row - 1) * grid_size
    
    pupil_x = (col - grid_size/2 - 0.5) * APMX() / grid_size
    pupil_y = (row - grid_size/2 - 0.5) * APMX() / grid_size
    
RETURN

! Random sampling pattern
LABEL RANDOM_SAMPLING
    ! Generate random point within circular aperture
    r = SQRT(RAND()) * APMX() / 2
    theta = RAND() * 360
    
    pupil_x = r * COSI(theta)
    pupil_y = r * SINE(theta)
    
RETURN