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
# Kokkos
set(Kokkos_DIR ${CMAKE_SOURCE_DIR}/includes/matar/install/kokkos/lib/cmake/kokkos)
# MATAR headers
set(matar_headers "${CMAKE_SOURCE_DIR}/includes/matar/install/matar/include")
# pprint headers
set(pprint_headers "${CMAKE_SOURCE_DIR}/includes/pprint/include")
# argparse header
set(argparse_headers "${CMAKE_SOURCE_DIR}/includes/argparse/include")

# ... ensure MPI is used
find_package(MPI REQUIRED)
find_package(VTK COMPONENTS CommonCore
                            CommonDataModel
                            IOXML
                 REQUIRED
)
find_package(Kokkos REQUIRED PATHS 
             ${CMAKE_SOURCE_DIR}/includes/matar/install/kokkos)
find_package(yaml-cpp REQUIRED PATHS
             ${CMAKE_SOURCE_DIR}/includes/yaml-cpp/install)

# ... Set compiler settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -debug -traceback")

file(GLOB sources "${src}/*.cpp") 
file(GLOB output "${src}/IO/*.cpp")
file(GLOB meshing "${src}/mesh/*.cpp")

list(APPEND cxx_sources ${sources}
                        ${output}
                        ${meshing}
)

# ... compile and create the executable
add_executable(tuskan ${cxx_sources})

# ... set compiler flag so that the correct MATAR headers are used
target_compile_definitions(tuskan PUBLIC HAVE_KOKKOS)

# ... include custom headers
target_include_directories(tuskan PUBLIC "${CMAKE_SOURCE_DIR}/includes")
# ... include MATAR headers
target_include_directories(tuskan PUBLIC ${matar_headers})
# ... include pprint headers
target_include_directories(tuskan PUBLIC ${pprint_headers})
# ... include argparse header
target_include_directories(tuskan PUBLIC ${argparse_headers})

# ... link the libraries
target_link_libraries(tuskan PRIVATE Kokkos::kokkos 
                                     yaml-cpp::yaml-cpp
                                     ${VTK_LIBRARIES}
                                     ${MPI_CXX_LIBRARIES}
)
# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS tuskan
  MODULES ${VTK_LIBRARIES}
)

# ... Install so that you can use anywhere
install(TARGETS tuskan DESTINATION bin)
