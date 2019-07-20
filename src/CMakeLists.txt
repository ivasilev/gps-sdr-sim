########################################################################
# Build libraries and executables
########################################################################
add_executable(gpssdrsim
	gpssim.h
	gpssim.c
)

if(MSVC)
	# needs CMake 3.1 but Windows builds should have that
	target_sources(gpssdrsim PRIVATE src/getopt.c)
endif()

if("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
    # silence this warning for now
    ADD_DEFINITIONS(-Wno-large-by-value-copy)
endif()



# Explicitly say that we want C99
set_property(TARGET gpssdrsim PROPERTY C_STANDARD 99)

