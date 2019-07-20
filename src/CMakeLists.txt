########################################################################
# Build libraries and executables
########################################################################
add_executable(rtl_433
	gpssim.h
	gpssim.c
)

if(MSVC)
	# needs CMake 3.1 but Windows builds should have that
	target_sources(rtl_433 PRIVATE getopt/getopt.c)
endif()

# Explicitly say that we want C99
set_property(TARGET rtl_433 PROPERTY C_STANDARD 99)

