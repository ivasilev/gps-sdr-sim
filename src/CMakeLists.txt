########################################################################
# Build libraries and executables
########################################################################
add_executable(gpssim
    constants.h
    generic_funcs.h
    generic_funcs.cpp
    gpstime.h
    gpstime.cpp
    ephemeris.h
    ephemeris.cpp
	gpssim.h
	gpssim.cpp
)

if(MSVC)
	# needs CMake 3.1 but Windows builds should have that
	target_sources(gpssim PRIVATE src/getopt.c)
endif()

if("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
    # silence this warning for now
    ADD_DEFINITIONS(-Wno-large-by-value-copy)
endif()



# Explicitly say that we want C99
set_property(TARGET gpssim PROPERTY CPP_STANDARD 14)

