#///////////////////////////////////////////////////////////////////////////
#//
#//    Copyright 2010
#//
#//    This file is part of rootpwa
#//
#//    rootpwa is free software: you can redistribute it and/or modify
#//    it under the terms of the GNU General Public License as published by
#//    the Free Software Foundation, either version 3 of the License, or
#//    (at your option) any later version.
#//
#//    rootpwa is distributed in the hope that it will be useful,
#//    but WITHOUT ANY WARRANTY; without even the implied warranty of
#//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#//    GNU General Public License for more details.
#//
#//    You should have received a copy of the GNU General Public License
#//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
#//
#///////////////////////////////////////////////////////////////////////////
#//-------------------------------------------------------------------------
#//
#// Description:
#//      cmake module for finding ROOT installation
#//      requires root-config to be in PATH
#//      based on AliRoots's FindROOT.cmake (r41015)
#//      in https://alisoft.cern.ch/AliRoot/trunk/cmake/modules
#//
#//      following variables are defined:
#//      ROOT_CONFIG_EXECUTABLE - path to root-config program
#//      ROOTSYS                - path to root installation directory
#//      ROOT_TARGET            - target architecture
#//      ROOT_F77               - Fortran complier used building ROOT
#//      ROOT_CC                - C complier used building ROOT
#//      ROOT_CPP               - C++ complier used building ROOT
#//      ROOT_VERSION           - ROOT version e.g. 5.34.0
#//      ROOT_MAJOR_VERSION     - ROOT major version e.g. 5
#//      ROOT_MINOR_VERSION     - ROOT minor version e.g. 34
#//      ROOT_PATCH_VERSION     - ROOT patch version e.g. 0
#//      ROOT_SVN_REVISION      - ROOT subversion revision
#//      ROOT_BIN_DIR           - ROOT executable directory
#//      ROOT_INCLUDE_DIR       - ROOT header directory
#//      ROOT_LIBRARY_DIR       - ROOT library directory
#//      ROOT_LIBRARIES         - linker flags for ROOT libraries
#//      ROOT_AUX_LIBRARIES     - linker flags for auxiliary libraries
#//      ROOT_EVE_LIBRARIES     - linker flags for EVE libraries
#//      ROOTCINT_EXECUTABLE    - path to rootcint program
#//      ROOT_LIBS              - list of ROOT library files
#//
#//      Example usage:
#//          find_package(ROOT 5.26 REQUIRED Minuit2)
#//
#//
#//      The module also provides a function to generate ROOT dictionaries.
#//      Example usage:
#//          set(ROOTPWA_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/someDict.cc)  # set dictionary path
#//          root_generate_dictionary(
#//            "${ROOTPWA_DICTIONARY}"            # path to dictionary to generate
#//            "${INCLUDE_DIR1};${INCLUDE_DIR2}"  # list of includes
#//            "class1.h;class2.h;class3.h"       # list of classes to process
#//            "someLinkDef.h"                    # ROOT linkDef file
#//          )
#//          set(SOURCES ${SOURCES} ${ROOTPWA_DICTIONARY})  # append dictionary to sources
#//
#//
#// Author List:
#//      Boris Grube          TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


set(ROOT_FOUND        FALSE)
set(ROOT_ERROR_REASON "")
set(ROOT_DEFINITIONS  "")
set(ROOT_LIBS)


find_program(ROOT_CONFIG_EXECUTABLE root-config)
if(NOT ROOT_CONFIG_EXECUTABLE)
	set(ROOT_ERROR_REASON "${ROOT_ERROR_REASON} Cannot find root-config executable in path. Make sure ROOT is setup correctly.")
else()

	set(ROOT_FOUND TRUE)

	execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --prefix
		OUTPUT_VARIABLE ROOTSYS
		OUTPUT_STRIP_TRAILING_WHITESPACE)

	execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --arch
		OUTPUT_VARIABLE ROOT_TARGET
		OUTPUT_STRIP_TRAILING_WHITESPACE)

	execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --f77
		OUTPUT_VARIABLE ROOT_F77
		OUTPUT_STRIP_TRAILING_WHITESPACE)

	execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --cc
		OUTPUT_VARIABLE ROOT_CC
		OUTPUT_STRIP_TRAILING_WHITESPACE)

	execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --cxx
		OUTPUT_VARIABLE ROOT_CPP
		OUTPUT_STRIP_TRAILING_WHITESPACE)

	execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --version
		OUTPUT_VARIABLE ROOT_VERSION
		OUTPUT_STRIP_TRAILING_WHITESPACE)

	execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --svn-revision
		OUTPUT_VARIABLE ROOT_SVN_REVISION
		OUTPUT_STRIP_TRAILING_WHITESPACE)

	execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --bindir
		OUTPUT_VARIABLE ROOT_BIN_DIR
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	if(NOT EXISTS "${ROOT_BIN_DIR}")
		set(ROOT_FOUND FALSE)
		set(ROOT_ERROR_REASON "${ROOT_ERROR_REASON} ROOT executable directory '${ROOT_BIN_DIR}' does not exist.")
	endif()

	execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --incdir
		OUTPUT_VARIABLE ROOT_INCLUDE_DIR
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	if(NOT EXISTS "${ROOT_INCLUDE_DIR}")
		set(ROOT_FOUND FALSE)
		set(ROOT_ERROR_REASON "${ROOT_ERROR_REASON} ROOT include directory '${ROOT_INCLUDE_DIR}' does not exist.")
	endif()

	execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --libdir
		OUTPUT_VARIABLE ROOT_LIBRARY_DIR
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	if(NOT EXISTS "${ROOT_LIBRARY_DIR}")
		set(ROOT_FOUND FALSE)
		set(ROOT_ERROR_REASON "${ROOT_ERROR_REASON} ROOT library directory '${ROOT_LIBRARY_DIR}' does not exist.")
	endif()

	execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --noauxlibs --glibs
		OUTPUT_VARIABLE ROOT_LIBRARIES
		OUTPUT_STRIP_TRAILING_WHITESPACE)

	# is this really needed? could one not handle this via the components?
  # in any case the EVE libs can be handled like the other ROOT libs
  # so this list should be merged with ROOT_LIBRARIES and duplicates removed
	# the extra loop over the EVE libs is superfluous
	execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --noauxlibs --evelibs
		OUTPUT_VARIABLE ROOT_EVE_LIBRARIES
		OUTPUT_STRIP_TRAILING_WHITESPACE)

	execute_process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --auxlibs
		OUTPUT_VARIABLE ROOT_AUX_LIBRARIES
		OUTPUT_STRIP_TRAILING_WHITESPACE)

	find_program(ROOTCINT_EXECUTABLE rootcint)
	if(NOT ROOTCINT_EXECUTABLE)
		set(ROOT_FOUND FALSE)
		set(ROOT_ERROR_REASON "${ROOT_ERROR_REASON} Cannot find rootcint. Make sure ROOT is setup correctly.")
	endif()

	# parse version string
	string(REGEX REPLACE "^([0-9]+)\\.[0-9][0-9]+\\/[0-9][0-9]+.*" "\\1"
		ROOT_MAJOR_VERSION "${ROOT_VERSION}")
	string(REGEX REPLACE "^[0-9]+\\.([0-9][0-9])+\\/[0-9][0-9]+.*" "\\1"
		ROOT_MINOR_VERSION "${ROOT_VERSION}")
	string(REGEX REPLACE "^[0-9]+\\.[0-9][0-9]+\\/([0-9][0-9]+).*" "\\1"
		ROOT_PATCH_VERSION "${ROOT_VERSION}")
	set(ROOT_VERSION "${ROOT_MAJOR_VERSION}.${ROOT_MINOR_VERSION}.${ROOT_PATCH_VERSION}")
	# compare version
	if(ROOT_FIND_VERSION_EXACT)
		if(NOT ROOT_VERSION VERSION_EQUAL ROOT_FIND_VERSION)
			set(ROOT_FOUND FALSE)
			set(ROOT_ERROR_REASON "${ROOT_ERROR_REASON} ROOT version ${ROOT_VERSION} does not match requested version ${ROOT_FIND_VERSION}.")
		endif()
	else()
		if(ROOT_VERSION VERSION_LESS ROOT_FIND_VERSION)
			set(ROOT_FOUND FALSE)
			set(ROOT_ERROR_REASON "${ROOT_ERROR_REASON} ROOT version ${ROOT_VERSION} is lower than requested version ${ROOT_FIND_VERSION}.")
		endif()
	endif()

endif()


# generate list of ROOT libraries
if(ROOT_FOUND)

	# create list of internal libraries from root-config output
	set(_LIBRARY_NAMES)
	set(_EXTERNAL_ZLIB)
	separate_arguments(ROOT_LIBRARIES)
	# remove first -L entry
	list(REMOVE_AT ROOT_LIBRARIES 0)
	# loop over -l entries
	foreach(_LIBRARY ${ROOT_LIBRARIES})
		# extract library name from compiler flag and append to list
		string(REGEX REPLACE "^-.(.*)$" "\\1" _LIBNAME "${_LIBRARY}")
		# workaround for root-config inconsistency: if ROOT is built with --disable-builtin-zlib
		# root-config returns the flag for the external zlib together with the internal libraries
		if(_LIBNAME STREQUAL "z")
			set(_EXTERNAL_ZLIB "-lz")
		else()
			list(APPEND _LIBRARY_NAMES ${_LIBNAME})
		endif()
		unset(_LIBRARY)
		unset(_LIBNAME)
	endforeach()

	# append components
	list(REMOVE_DUPLICATES ROOT_FIND_COMPONENTS)
	if(ROOT_FIND_COMPONENTS)
		set(_LIBRARY_NAMES "${_LIBRARY_NAMES};${ROOT_FIND_COMPONENTS}")
	endif()

	# check whether libraries exist
	foreach(_LIBNAME ${_LIBRARY_NAMES})
		find_library(_ROOT_LIB_${_LIBNAME}
			NAMES ${_LIBNAME}
			PATHS ${ROOT_LIBRARY_DIR}
			NO_DEFAULT_PATH)
		if(NOT _ROOT_LIB_${_LIBNAME})
			set(ROOT_FOUND FALSE)
			set(ROOT_ERROR_REASON "${ROOT_ERROR_REASON} Cannot find ROOT library '${_LIBNAME}' in '${ROOT_LIBRARY_DIR}'.")
		else()
			list(APPEND ROOT_LIBS ${_ROOT_LIB_${_LIBNAME}})
		endif()
		unset(_ROOT_LIB_${_LIBNAME})
	endforeach()
	unset(_LIBNAME)

	# create list of external libraries from root-config output
	separate_arguments(ROOT_AUX_LIBRARIES)
	# append external zlib to auxiliary libraries
	if(_EXTERNAL_ZLIB)
		list(APPEND ROOT_AUX_LIBRARIES ${_EXTERNAL_ZLIB})
	endif()
	unset(_EXTERNAL_ZLIB)
	# loop over -l entries
	foreach(_LIBRARY ${ROOT_AUX_LIBRARIES})
		# extract library name from compiler flag
		string(REGEX MATCH "^-l(.*)$" _LIBNAME "${_LIBRARY}")
		if(_LIBNAME)
			string(REGEX REPLACE "^-.(.*)$" "\\1" _LIBNAME "${_LIBNAME}")
			# check whether libraries exist
			find_library(_AUX_LIB_${_LIBNAME}
				NAMES ${_LIBNAME})
			if(NOT _AUX_LIB_${_LIBNAME})
				set(ROOT_FOUND FALSE)
				set(ROOT_ERROR_REASON "${ROOT_ERROR_REASON} Cannot find ROOT library '${_LIBNAME}'.")
			else()
				list(APPEND ROOT_LIBS ${_AUX_LIB_${_LIBNAME}})
			endif()
			unset(_AUX_LIB_${_LIBNAME})
		endif()
		unset(_LIBNAME)
	endforeach()
	unset(_LIBRARY)

	# create list of external libraries from root-config output, for evelibs this time
	separate_arguments(ROOT_EVE_LIBRARIES)
	# remove first -L entry
	list(REMOVE_AT ROOT_EVE_LIBRARIES 0)
	# loop over -l entries
	foreach(_LIBRARY ${ROOT_EVE_LIBRARIES})
		# extract library name from compiler flag
		string(REGEX MATCH "^-l...(.*)$" _LIBNAME "${_LIBRARY}")
		if(_LIBNAME)
			string(REGEX REPLACE "^-.(.*)$" "\\1" _LIBNAME "${_LIBNAME}")
			# check whether libraries exist
			find_library(_EVE_LIB_${_LIBNAME}
				NAMES ${_LIBNAME}
				HINTS ${ROOT_LIBRARY_DIR})
			if(NOT _EVE_LIB_${_LIBNAME})
				set(ROOT_FOUND FALSE)
				set(ROOT_ERROR_REASON "${ROOT_ERROR_REASON} Cannot find ROOT library '${_LIBNAME}'.")
			else()
				list(APPEND ROOT_LIBS ${_EVE_LIB_${_LIBNAME}})
				list(APPEND ROOT_LIBRARIES -l${_LIBNAME})
			endif()
			unset(_EVE_LIB_${_LIBNAME})
		endif()
		unset(_LIBNAME)
	endforeach()

endif()


# make variables changeable
mark_as_advanced(
	ROOT_INCLUDE_DIR
	ROOT_LIBRARY_DIR
	ROOT_LIBRARIES
	ROOT_LIBS
	ROOT_DEFINITIONS
	)


# report result
if(ROOT_FOUND)
	message(STATUS "Found ROOT version ${ROOT_VERSION} r${ROOT_SVN_REVISION} in '${ROOTSYS}'.")
	message(STATUS "Using ROOT include directory '${ROOT_INCLUDE_DIR}'.")
	message(STATUS "Using ROOT library directory '${ROOT_LIBRARY_DIR}'.")
	message(STATUS "Using ROOT libraries ${ROOT_LIBRARIES}.")
	message(STATUS "Using ROOT additional components '${ROOT_FIND_COMPONENTS}'.")
else()
	if(ROOT_FIND_REQUIRED)
		message(FATAL_ERROR "Unable to find requested ROOT installation:${ROOT_ERROR_REASON}")
	else()
		if(NOT ROOT_FIND_QUIETLY)
			message(STATUS "ROOT version ${ROOT_FIND_VERSION}+ was not found:${ROOT_ERROR_REASON}")
		endif()
	endif()
endif()


# function that generates ROOT dictionary
function(root_generate_dictionary DICT_FILE INCLUDE_DIRS HEADER_FILES LINKDEF_FILE)

	if(DEBUG_OUTPUT)
		message(STATUS "root_generate_dictionary was called with the following arguments:
        DICT_FILE    = '${DICT_FILE}'
        INCLUDE_DIRS = '${INCLUDE_DIRS}'
        HEADER_FILES = '${HEADER_FILES}'
        LINKDEF_FILE = '${LINKDEF_FILE}'")
	endif()

	if(NOT ROOT_FOUND)
		message(FATAL_ERROR "Impossible to generate dictionary '${DICT_FILE}', "
			"because no ROOT installation was found.")
	endif()

	# prepare command line argument for compiler definitions (put -D in front)
	set(_DEFINITIONS)
	get_property(_DEFS DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY COMPILE_DEFINITIONS)
	foreach(_DEF ${_DEFS})
		set(_DEFINITIONS "${_DEFINITIONS} -D${_DEF}")
	endforeach()
	unset(_DEF)
	separate_arguments(_DEFINITIONS)

	# prepare command line argument for include directories (put -I in front)
	set(_INCLUDES)
	foreach(_FILE ${INCLUDE_DIRS})
		set(_INCLUDES ${_INCLUDES} -I${_FILE})
	endforeach()
	unset(_FILE)

	# strip paths from header file names
	set(_HEADERS)
	foreach(_FILE ${HEADER_FILES})
		get_filename_component(_NAME ${_FILE} NAME)
		set(_HEADERS ${_HEADERS} ${_NAME})
	endforeach()
	unset(_FILE)

	# add dictionary header file to output files
	string(REGEX REPLACE "^(.*)\\.(.*)$" "\\1.h" _DICT_HEADER "${DICT_FILE}")
	set(OUTPUT_FILES ${DICT_FILE} ${_DICT_HEADER})
	if(DEBUG_OUTPUT)
		message(STATUS "root_generate_dictionary will create output files '${OUTPUT_FILES}'")
	endif()

	add_custom_command(OUTPUT ${OUTPUT_FILES}
		COMMAND ${ROOTCINT_EXECUTABLE}
		ARGS -f ${DICT_FILE} -c -DHAVE_CONFIG_H ${_DEFINITIONS} ${_INCLUDES} ${_HEADERS} ${LINKDEF_FILE}
		DEPENDS ${HEADER_FILES} ${LINKDEF_FILE}
		)
	if(DEBUG_OUTPUT)
		message(STATUS "root_generate_dictionary will execute
        '${ROOTCINT_EXECUTABLE} -f ${DICT_FILE} -c -DHAVE_CONFIG_H ${_DEFINITIONS} ${_INCLUDES} ${_HEADERS} ${LINKDEF_FILE}'")
	endif()
	unset(_DEFINITIONS)
	unset(_INCLUDES)
	unset(_HEADERS)

endfunction(root_generate_dictionary)
