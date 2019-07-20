########################################################################
# Project setup
########################################################################
cmake_minimum_required(VERSION 2.6)
# Fix behavior of CMAKE_C_STANDARD when targeting macOS.
if(POLICY CMP0025)
    cmake_policy(SET CMP0025 NEW)
endif()

project(gpsSdrSim)

########################################################################
# Compiler specific setup
########################################################################
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 99)
if(("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_C_COMPILER_ID}" MATCHES "Clang") AND NOT WIN32)
    ADD_DEFINITIONS(-Wall)
    ADD_DEFINITIONS(-Wextra)
    ADD_DEFINITIONS(-Wsign-compare)
    ADD_DEFINITIONS(-Wno-unused)
    ADD_DEFINITIONS(-Wno-unused-parameter)
    ADD_DEFINITIONS(-Wno-missing-field-initializers)
    ADD_DEFINITIONS(-std=c99)
    ADD_DEFINITIONS(-pedantic)
    # for strdup, setenv, use either
    #ADD_DEFINITIONS(-D_POSIX_C_SOURCE=200809) # does not work with uClibc
    ADD_DEFINITIONS(-D_GNU_SOURCE)
    #http://gcc.gnu.org/wiki/Visibility
    add_definitions(-fvisibility=hidden)

    # CMake Release default for GCC/Clang is "-O3 -DNDEBUG"
    # set(CMAKE_C_FLAGS_RELEASE -O2)
    # CMake Debug default for GCC/Clang is "-g -DNDEBUG"
    # set(CMAKE_C_FLAGS_DEBUG -g3 -O0)
    # make use of ASAN
    set(CMAKE_C_FLAGS_DEBUG "-ggdb -fsanitize=undefined -fsanitize=address -fno-omit-frame-pointer")
endif()
if("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
    # make sure we don't accidentally copy more than an int
    ADD_DEFINITIONS(-Wlarge-by-value-copy=8)
endif()

# Shut MSVC up about strdup and strtok
if(MSVC)
    ADD_DEFINITIONS(-D_CRT_NONSTDC_NO_DEPRECATE)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
    ADD_DEFINITIONS(-DNOMINMAX)
endif()

# Fix printf %zu
if(MINGW)
    add_definitions(-D__USE_MINGW_ANSI_STDIO)
endif()

########################################################################
# Find build dependencies
########################################################################
# cmake -DCMAKE_BUILD_TYPE=Profile ..
# CPUPROFILE=prof.out ./src/rtl_433 ...
# pprof -text ./src/rtl_433 prof.out
if("${CMAKE_BUILD_TYPE}" STREQUAL "Profile")
    message(STATUS "Build type set to Profile. Linking GPerfTools.")
    find_package(Gperftools REQUIRED)
    include_directories(${GPERFTOOLS_INCLUDE_DIR})
    list(APPEND SDR_LIBRARIES ${GPERFTOOLS_LIBRARIES} -Wl,-no_pie)
    ADD_DEFINITIONS(-g)
    ADD_DEFINITIONS(-fno-builtin-malloc)
    ADD_DEFINITIONS(-fno-builtin-calloc)
    ADD_DEFINITIONS(-fno-builtin-realloc)
    ADD_DEFINITIONS(-fno-builtin-free)
endif()

########################################################################
# Setup the include and linker paths
########################################################################
include_directories(
    BEFORE
    ${CMAKE_SOURCE_DIR}/include
)

########################################################################
# Add subdirectories
########################################################################
add_subdirectory(src)


