cmake_minimum_required(VERSION 3.15.0)

project(vmec LANGUAGES Fortran)

project (VMEC Fortran)

set(CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "ON")

# Set a directories to build all binary files.
set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)
set (LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib)

# Fortran specific settings. The first setting tells the compiler to use the C preprocessor.
# The second places a common directory to place all of the module files.
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

find_package (BLAS REQUIRED)

set (vmec_sources "")
add_subdirectory(src)

add_library (vmec STATIC ${vmec_sources})
set_target_properties (vmec PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/vmec)
target_include_directories(vmec PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY}/vmec)

if (BLAS_FOUND)
	target_link_libraries (vmec PUBLIC ${BLAS_LIBRARIES})
endif ()

get_target_property(VMEC_LINK_LIB vmec LINK_LIBRARIES)
message(STATUS "vmec linked libraries are ${VMEC_LINK_LIB}")

# Define an executable and link all libraries.
add_executable (xvmec ${CMAKE_CURRENT_SOURCE_DIR}/src/vmec.f90)
target_link_libraries (xvmec vmec)
