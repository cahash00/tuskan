# CMakeLists.txt in src/
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "*.cpp")
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS "*.h")

# Add the include directories for all subdirectories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Add sources to the target
add_library(BCs ${SOURCES} ${HEADERS})
target_include_directories(BCs PRIVATE ${matar_headers})
target_include_directories(BCs PRIVATE ${src}/common)
target_include_directories(BCs PRIVATE ${src}/IO)
target_include_directories(BCs PRIVATE ${YAML_CPP_INCLUDE_DIR})
