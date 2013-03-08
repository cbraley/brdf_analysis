cmake_minimum_required(VERSION 2.8)

# High Level setup ------------------------------------------------------------
project(BRDF_Analyze) # Project name.
set(BINARY_NAME bin/brdf_analyze) # Output location of the final C++ binary.
set(CMAKE_BUILD_TYPE Release) # Choose from {Debug, Release}.  Affects
# optimization & debug symbols.
set(CMAKE_MODULE_PATH CMakeModules) # Location of the Find.*.cmake files.

# Find source and header files ------------------------------------------------

# Lists of source files.
file(GLOB SRCS "src/*cpp")

# Print debugging info --------------------------------------------------------
message(STATUS "Source files:")
foreach(FILE ${SRCS})
    message(STATUS "\t${FILE}")
endforeach()

# Add the binary --------------------------------------------------------------
add_executable(${BINARY_NAME}  ${SRCS})

# GCC flags -------------------------------------------------------------------
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -W -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g0 -msse2")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3")

# pthreads library ------------------------------------------------------------
# TODO(Colin): pthreads is currently not used...

set(CMAKE_THREAD_PREFER_PTHREAD 1) #In the presence of choices, use pthreads.
find_package(Threads REQUIRED)
if(${Threads_FOUND})
    message("Successfully found threading library.")
    target_link_libraries(${BINARY_NAME} ${CMAKE_THREAD_LIBS_INIT})
else()
    message(FATAL_ERROR "Could not find a suitable threads library(usually pthreads)!")
endif()

# BLAS  ------------------------------------------------------------------------
# Note - BLAS is not needed if you define DONT_USE_BLAS
find_package(BLAS REQUIRED)
if(${BLAS_FOUND})
    message("Successfully found blas with vendor ${BLA_VENDOR}.")
    target_link_libraries(${BINARY_NAME} ${BLAS_LIBRARIES} ${BLAS95_LIBRARIES})
else()
    message(FATAL_ERROR "Could not find BLAS; try defining DONT_USE_BLAS.")
endif()


