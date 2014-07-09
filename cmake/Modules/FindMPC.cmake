# - Try to find the MPC libraries
# This module defines:
# MPC_FOUND - system has MPC lib
# MPC_INCLUDE_DIR - the MPC include directory
# MPC_LIBRARIES_DIR - directory where the MPC libraries are located
# MPC_LIBRARIES - Link these to use MPC
# MPC_IN_CGAL_AUXILIARY - TRUE if the MPC found is the one distributed with CGAL in the auxiliary folder

# TODO: support MacOSX

include(FindPackageHandleStandardArgs)

if(MPC_INCLUDE_DIR)
  set(MPC_in_cache TRUE)
else()
  set(MPC_in_cache FALSE)
endif()
if(NOT MPC_LIBRARIES)
  set(MPC_in_cache FALSE)
endif()

# Is it already configured?
if (MPC_in_cache)

  set(MPC_FOUND TRUE)

else()

  find_path(MPC_INCLUDE_DIR
            NAMES mpc.h
            HINTS ENV MPC_INC_DIR
                  ENV MPC_DIR
                  ENV MPC_ROOT
            PATHS ENV MPC_ROOT ENV MPC_DIR ENV MPC_INC_DIR /usr/include /usr/local/include ENV MPC_ROOT
            PATH_SUFFIXES include
   DOC "The directory containing the MPC header files"
           )

  find_library(MPC_LIBRARIES NAMES mpc libmpc
    HINTS ENV MPC_LIB_DIR
          ENV MPC_DIR
          ENV MPC_ROOT
    PATHS ENV MPC_ROOT ENV MPC_DIR ENV MPC_LIB_DIR /usr/lib64 /usr/local/lib64 /usr/lib /usr/lib64
    PATH_SUFFIXES lib lib64
    DOC "Path to the MPC library"
    )

  if ( MPC_LIBRARIES )
    get_filename_component(MPC_LIBRARIES_DIR ${MPC_LIBRARIES} PATH CACHE )
  endif()

  # Attempt to load a user-defined configuration for MPC if couldn't be found
  if ( NOT MPC_INCLUDE_DIR OR NOT MPC_LIBRARIES_DIR )
    include( MPCConfig OPTIONAL )
  endif()

  find_package_handle_standard_args(MPC "DEFAULT_MSG" MPC_LIBRARIES MPC_INCLUDE_DIR)

endif()