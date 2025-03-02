# ... set minimum CMake version required
cmake_minimum_required(VERSION 3.16.0)

# ... set project details
project(tuskan VERSION 0.0
   DESCRIPTION "TUSKAN Two-Phase Flow Solver")

# ... Enable Fortran and CXX for the project
enable_language(CXX)

# ... set CMake Build Type
set(CMAKE_BUILD_TYPE Release)
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Develop" "Debug")
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries")

# ... set final install directory
set(CMAKE_INSTALL_PREFIX
    "${CMAKE_BINARY_DIR}/install"
    CACHE STRING
    "Final install location." FORCE
)

# ... Set paths and directories
# source files
set(src ${CMAKE_SOURCE_DIR}/src)
# include path
set(includes ${CMAKE_SOURCE_DIR}/includes)
# Kokkos
set(Kokkos_DIR ${includes}/matar/install/kokkos/lib/cmake/kokkos)
# MATAR headers
set(matar_headers "${includes}/matar/install/matar/include")
# pprint headers
set(pprint_headers "${includes}/pprint/include")
# argparse header
set(argparse_headers "${includes}/argparse/include")
# spdlog headers
set(spdlog_headers "${includes}/spdlog/include")

# ... ensure MPI is used
find_package(MPI REQUIRED)
find_package(Kokkos REQUIRED PATHS 
             ${CMAKE_SOURCE_DIR}/includes/matar/install/kokkos)
find_package(yaml-cpp REQUIRED PATHS
             ${CMAKE_SOURCE_DIR}/includes/yaml-cpp/build/install)

# ... Set compiler settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -debug -traceback")

file(GLOB sources "${src}/*.cpp") 

list(APPEND cxx_sources ${sources})

# ... add the custom libraries
add_subdirectory(${src}/mesh)
add_subdirectory(${src}/common)
add_subdirectory(${src}/IO)
add_subdirectory(${src}/utils)
add_subdirectory(${src}/solver)

# ... compile and create the executable
add_executable(tuskan ${cxx_sources})

# ... set compiler flag so that the correct MATAR headers are used
target_compile_definitions(tuskan PUBLIC HAVE_KOKKOS)

# ... include custom headers
target_include_directories(tuskan PRIVATE "${src}/mesh")
target_include_directories(tuskan PRIVATE "${src}/common")
target_include_directories(tuskan PRIVATE "${src}/IO")
target_include_directories(tuskan PRIVATE "${src}/utils")
target_include_directories(tuskan PRIVATE "${src}/solver")
# ... include MATAR headers
target_include_directories(tuskan PRIVATE ${matar_headers})
# ... include pprint headers
target_include_directories(tuskan PRIVATE ${pprint_headers})
# ... include argparse header
target_include_directories(tuskan PRIVATE ${argparse_headers})
# ... include the logging headers
target_include_directories(tuskan PRIVATE ${spdlog_headers})

# ... link the libraries
target_link_libraries(tuskan PRIVATE Kokkos::kokkos 
                                     yaml-cpp::yaml-cpp
                                     ${VTK_LIBRARIES}
                                     ${MPI_CXX_LIBRARIES}
                                     meshLib
                                     IOLib
                                     paramsLib
                                     generalUtils
                                     solvers
)

# ... Install so that you can use anywhere
install(TARGETS tuskan DESTINATION bin)
