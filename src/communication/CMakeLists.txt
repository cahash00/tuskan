# CMakeLists.txt in src/
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "*.cpp")
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS "*.h")

# Add the include directories for all subdirectories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Add sources to the target
add_library(comm ${SOURCES} ${HEADERS})
target_include_directories(comm PRIVATE ${pprint_headers})
target_include_directories(comm PRIVATE ${matar_headers})
target_include_directories(comm PRIVATE ${src}/common)
target_link_libraries(comm PRIVATE MPI::MPI_C)
