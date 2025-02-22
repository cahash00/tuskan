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
# Kokkos
set(src ${CMAKE_SOURCE_DIR}/src)
set(Kokkos_DIR ${CMAKE_SOURCE_DIR}/includes/matar/install/kokkos/lib/cmake/kokkos)
# MATAR headers
set(matar_headers "${CMAKE_SOURCE_DIR}/includes/matar/install/matar/include")
set(pprint_headers "${CMAKE_SOURCE_DIR}/includes/pprint/include")

# ... ensure MPI is used
find_package(MPI REQUIRED)
find_package(Kokkos REQUIRED PATHS 
             ${CMAKE_SOURCE_DIR}/includes/matar/install/kokkos)

# ... Set compiler settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -debug -traceback")

file(GLOB cxx_sources "${src}/*.cpp") 

# ... compile and create the executable
add_executable(tuskan ${cxx_sources})

# ... set compiler flag so that the correct MATAR headers are used
target_compile_definitions(tuskan PUBLIC HAVE_KOKKOS)

# ... include MATAR headers
target_include_directories(tuskan PUBLIC ${matar_headers})
# ... include pprint headers
target_include_directories(tuskan PUBLIC ${pprint_headers})

# ... link the libraries
target_link_libraries(tuskan PRIVATE Kokkos::kokkos
   ${MPI_CXX_LIBRARIES}
)

# ... Install so that you can use anywhere
install(TARGETS tuskan DESTINATION bin)
