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
 "${PROJECT_SOURCE_DIR}/src/odexpConfig.h.in"
 "${PROJECT_BINARY_DIR}/odexpConfig.h"
 )

# ============================================================================
# add the binary tree to the search path for include files so that we will find odexpConfig.h
# ============================================================================
include_directories("${PROJECT_BINARY_DIR}")

# ============================================================================
# Find GSL
# ============================================================================
find_package(GSL REQUIRED)

# ============================================================================
# add library target 
# ============================================================================
add_library(odexp SHARED
  src/datastruct.c
  src/main.c
  src/methods_odexp.c
  src/rand_gen.c
  src/utils_odexp.c
)
target_link_libraries(odexp GSL::gsl GSL::gslcblas)
target_link_libraries(odexp /usr/local/opt/readline/lib/libreadline.dylib)
target_link_libraries(odexp /usr/local/opt/readline/lib/libhistory.dylib)
set_target_properties(odexp PROPERTIES VERSION ${PROJECT_VERSION})

# ============================================================================
# Declare public API of your library
# ============================================================================
set_target_properties(odexp PROPERTIES PUBLIC_HEADER src/odexp.h)

# ============================================================================
# private headers should not been installed
# ============================================================================
target_include_directories(odexp PRIVATE src)

# ============================================================================
# Get GNU standard installation directories (GNUInstallDirs module)
# ============================================================================
include(GNUInstallDirs)

# ============================================================================
# Installation dirs 
# ============================================================================

IF ( NOT DEFINED MANDIR )
  set( MANDIR ${CMAKE_INSTALL_PREFIX}/share/man/man3 )
ENDIF ( NOT DEFINED MANDIR )
message( "** Manual files will be installed in ${MANDIR}" )

IF ( NOT DEFINED BINDIR )
  set( BINDIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}  )
ENDIF ( NOT DEFINED BINDIR )
message( "** Executable files will be installed in ${BINDIR}" )

# ============================================================================
# Set compilation flags (general and then build-type specific)
# ============================================================================
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -Wextra -Wno-sign-compare -O0")
set(CMAKE_BUILD_TYPE Debug)


# ============================================================================
# declare files to install
# ============================================================================
install(TARGETS odexp
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# ============================================================================
# add a rule to expand @ macros (@ONLY ask to cmake to not expand variables 
# of the form ${VAR})
# ============================================================================
configure_file(src/odexp.pc.in odexp.pc @ONLY)

# ============================================================================
# Add custom command for manpage
# ============================================================================
ADD_CUSTOM_TARGET( manpages ALL )
ADD_CUSTOM_COMMAND( TARGET manpages
  COMMAND cp ${PROJECT_SOURCE_DIR}/src/odexp.man ${CMAKE_BINARY_DIR}/odexp.3
  OUTPUTS ${CMAKE_BINARY_DIR}/odexp.3
  DEPENDS ${CMAKE_BINARY_DIR}/odexp.3
)

# ============================================================================
# install generated files
# ============================================================================
install(FILES ${CMAKE_BINARY_DIR}/odexp.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)

install(FILES ${CMAKE_BINARY_DIR}/odexp.3 DESTINATION ${MANDIR}) 

install(PROGRAMS ${CMAKE_SOURCE_DIR}/odexp DESTINATION ${BINDIR} )

