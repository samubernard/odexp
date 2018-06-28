# ============================================================================
# Require minimal version of cmake
# ============================================================================
cmake_minimum_required (VERSION 3.6)

# ============================================================================
# Set project name and languages
# ============================================================================
project (odexp
        VERSION 1.0.0
        LANGUAGES C)

# ============================================================================
# Set project version number
# ============================================================================
set (odexp_VERSION_MAJOR 1)
set (odexp_VERSION_MINOR 0)

 
# ============================================================================
# configure a header file to pass some of the CMake settings to the source code
# ============================================================================
configure_file (
  "${PROJECT_SOURCE_DIR}/odexpConfig.h.in"
  "${PROJECT_BINARY_DIR}/odexpConfig.h"
  )

# ============================================================================
# add the binary tree to the search path for include files so that we will find TutorialConfig.h
# ============================================================================
include_directories("${PROJECT_BINARY_DIR}")

# ============================================================================
# add library target 
# ============================================================================
add_library(odexp SHARED
  src/dlist.c
  src/methods_odexp.c
  src/odexp.c
  src/rand_gen.c
  src/utils_odexp.c
)
set_target_properties(odexp PROPERTIES VERSION ${PROJECT_VERSION})

# ============================================================================
# Declare public API of your library
# ============================================================================
set_target_properties(odexp PROPERTIES PUBLIC_HEADER include/odexp.h)

target_include_directories(odexp PRIVATE src)







