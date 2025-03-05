# CMakeLists.txt in src/
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "*.cpp")
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS "*.h")

# Add the include directories for all subdirectories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Add sources to the target
add_library(IOLib ${SOURCES} ${HEADERS})

target_include_directories(IOLib PRIVATE ${argparse_headers})
target_include_directories(IOLib PRIVATE ${pprint_headers})
target_include_directories(IOLib PRIVATE ${matar_headers})
target_include_directories(IOLib PRIVATE ${src}/common)
# just include the YAML headers, you do not need to link again!
target_include_directories(IOLib PRIVATE ${YAML_CPP_INCLUDE_DIR})
