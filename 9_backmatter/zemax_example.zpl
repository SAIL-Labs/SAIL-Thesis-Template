! Simple ZEMAX ZPL Macro for optical analysis
! Basic system information and ray tracing

! Get system parameters
num_surfaces = NSUR()
num_wavelengths = NWAV()
num_fields = NFLD()

! Print system info
PRINT "System Analysis"
PRINT "==============="
PRINT "Surfaces: ", num_surfaces
PRINT "Wavelengths: ", num_wavelengths
PRINT "Fields: ", num_fields
PRINT ""

! Calculate basic parameters
focal_length = 1/PWAV(1)
aperture_diameter = APMX()
f_number = focal_length / aperture_diameter

PRINT "Focal Length: ", focal_length, " mm"
PRINT "Aperture: ", aperture_diameter, " mm" 
PRINT "F-number: f/", f_number
PRINT ""

! Simple ray tracing for each field
PRINT "Ray Tracing Results:"
PRINT "-------------------"

FOR field_num, 1, num_fields, 1
    field_x = FLDX(field_num)
    field_y = FLDY(field_num)
    
    ! Trace chief ray
    ray_status = RAYT(0, 0, field_x, field_y, 1)
    
    IF ray_status = 1
        image_x = RAYX(num_surfaces)
        image_y = RAYY(num_surfaces)
        image_z = RAYZ(num_surfaces)
        
        PRINT "Field ", field_num, ": ("
        PRINT image_x, ", ", image_y, ", ", image_z, ")"
    ELSE
        PRINT "Field ", field_num, ": Ray failed"
    ENDIF
NEXT field_num

PRINT ""
PRINT "Analysis complete"