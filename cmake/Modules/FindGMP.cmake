# - Try to find the GMP libraries
# This module defines:
# GMP_FOUND - system has GMP lib
# GMP_INCLUDE_DIR - the GMP include directory
# GMP_LIBRARIES_DIR - directory where the GMP libraries are located
# GMP_LIBRARIES - Link these to use GMP
# GMP_IN_CGAL_AUXILIARY - TRUE if the GMP found is the one distributed with CGAL in the auxiliary folder

# TODO: support MacOSX

include(FindPackageHandleStandardArgs)

if(GMP_INCLUDE_DIR)
  set(GMP_in_cache TRUE)
else()
  set(GMP_in_cache FALSE)
endif()
if(NOT GMP_LIBRARIES)
  set(GMP_in_cache FALSE)
endif()

# Is it already configured?
if (GMP_in_cache)

  set(GMP_FOUND TRUE)

else()

  find_path(GMP_INCLUDE_DIR
            NAMES gmp.h
            HINTS ENV GMP_INC_DIR
                  ENV GMP_DIR
                  ENV GMP_ROOT
            PATHS ENV GMP_ROOT ENV GMP_DIR ENV GMP_INC_DIR /usr/include /usr/local/include ENV GMP_ROOT
            PATH_SUFFIXES include
   DOC "The directory containing the GMP header files"
           )

  if (GMP_INCLUDE_DIR)
    # extract version
    file(READ "${GMP_INCLUDE_DIR}/gmp.h" _GMP_VERSION_FILE)

    string(REGEX REPLACE ".*#define __GNU_MP_VERSION[ \t]*([0-9]+)[ \t\r\n]*.*" "\\1" GMP_VERSION_MAJOR "${_GMP_VERSION_FILE}")
    string(REGEX REPLACE ".*#define __GNU_MP_VERSION_MINOR[ \t]*([0-9]+)[ \t\r\n]*.*" "\\1" GMP_VERSION_MINOR "${_GMP_VERSION_FILE}")
    string(REGEX REPLACE ".*#define __GNU_MP_VERSION_PATCHLEVEL[ \t]*([0-9]+)[ \t\r\n]*.*" "\\1" GMP_VERSION_PATCH "${_GMP_VERSION_FILE}")
    unset(_GMP_VERSION_FILE)
    set(GMP_VERSION "${GMP_VERSION_MAJOR}.${GMP_VERSION_MINOR}.${GMP_VERSION_PATCH}")
  endif()

  find_library(GMP_LIBRARIES NAMES gmp libgmp-10 libgmp.10 gmp.10 gmp-10
    HINTS ENV GMP_LIB_DIR
          ENV GMP_DIR
          ENV GMP_ROOT
    PATH_SUFFIXES lib lib64
    PATHS ENV GMP_ROOT ENV GMP_DIR ENV GMP_LIB_DIR /usr/lib64 /usr/local/lib64 /usr/lib /usr/lib64
    DOC "Path to the GMP library"
    )

  if ( GMP_LIBRARIES )
    get_filename_component(GMP_LIBRARIES_DIR ${GMP_LIBRARIES} PATH CACHE )
  endif()

  # Attempt to load a user-defined configuration for GMP if couldn't be found
  if ( NOT GMP_INCLUDE_DIR OR NOT GMP_LIBRARIES_DIR )
    include( GMPConfig OPTIONAL )
  endif()

  find_package_handle_standard_args(GMP "DEFAULT_MSG" GMP_LIBRARIES GMP_INCLUDE_DIR)

endif()

if(GMP_FOUND AND NOT TARGET GMP::GMP)
  add_library(GMP::GMP INTERFACE IMPORTED)
  set_target_properties(GMP::GMP PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}")
  target_link_libraries(GMP::GMP INTERFACE "${GMP_LIBRARIES}")
endif()
