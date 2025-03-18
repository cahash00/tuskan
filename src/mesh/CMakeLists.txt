# CMakeLists.txt in src/
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "*.cpp")
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS "*.h")

# Add the include directories for all subdirectories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Add sources to the target
add_library(meshLib ${SOURCES} ${HEADERS})
target_include_directories(meshLib PRIVATE ${matar_headers})
target_include_directories(meshLib PRIVATE ${pprint_headers})
target_include_directories(meshLib PRIVATE ${src}/common)
