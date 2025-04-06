# CMakeLists.txt in src/
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "*.cpp")
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS "*.h")

# Add the include directories for all subdirectories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Add sources to the target
add_library(solvers ${SOURCES} ${HEADERS})
target_include_directories(solvers PRIVATE ${matar_headers})
target_include_directories(solvers PRIVATE ${pprint_headers})
target_include_directories(solvers PRIVATE ${src}/common)
target_include_directories(solvers PRIVATE ${src}/BCs)
target_include_directories(solvers PRIVATE ${src}/IO)
target_include_directories(solvers PRIVATE ${YAML_CPP_INCLUDE_DIR})
