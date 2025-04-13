# CMakeLists.txt in src/
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "*.cpp")
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS "*.h")

# Add the include directories for all subdirectories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Add sources to the target
add_library(paramsLib ${SOURCES} ${HEADERS})  # Change to add_executable if needed
target_include_directories(paramsLib PRIVATE "${src}/communication")
target_link_libraries(paramsLib PRIVATE MPI::MPI_C)
